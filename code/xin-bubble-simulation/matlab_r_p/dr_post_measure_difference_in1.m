clc;
clear;
close all;
% Get screen size
% screenSize = get(groot,'ScreenSize');
% init_idx=10;
init_idx=[0,50,100,150];
% for i=1:numel(init_idx)
%     disp(init_idx(i))
% end
% % Set figure position
% figure('Position',[1 1 screenSize(3) screenSize(4)]);
% %%%%%%%%%%% MY-LBM %%%%%%%%%%%%
cd("..\")
% Define file names
% file_names = ["testcase2\pressure_100.00_101_0.31_VSM.dat", "testcase1\pressure_100.00_101_0.34_VSM.dat", "testcase1\pressure_400.00_401_0.31_VSM.dat", "testcase1\pressure_400.00_401_0.34_VSM.dat"];
parameters = jsondecode(fileread('input_matlab.json'));
radius=parameters.radius;
% radius(1)
total_time=parameters.mstep;
% file_names = cell(1, numel(parameters.mx)*numel(parameters.rho_l));
k=1;
file_names=cell(1,2);
legend_entries1=cell(1,numel(file_names));
% myidx = 1,2,3,4;
rad_lbm = {};
t_rad_lbm = {};
R_={};
t_={};
parameters.mx=1000;
parameters.lx=1001;
parameters.radius=30;
for i = 1:numel(parameters.mx)
    for j = 1:numel(parameters.rho_l)
%         for l = 2:numel(parameters.radius)
            file_names{k} = sprintf('testcase2\\pressure_%.0f_%d_%.2f_%d_VSM.dat', parameters.mx(i),parameters.lx(i), parameters.rho_l(j),parameters.radius);
%             legend_entries1{k} = sprintf('LBM-%d-%.2f-%d', parameters.lx(i), parameters.rho_l(j),parameters.radius(l));
%             legend_entries1{k} = sprintf('LBM-%.2f',parameters.rho_l(j));
            legend_entries1{k} = sprintf('R_{LBM}');
            k=k+1;
%         end
    end
end
% Loop over file names and plot data
for i = 1:numel(file_names)
    % Load data
%     if(i==myidx)
    data = load(file_names{i});

    % Calculate radii
    eval(sprintf('rad_lbm%d = 1./data(:,3);', i));
    rad_lbm{i} = 1./data(:,3);
    t_rad_lbm{i} = linspace(0,length(rad_lbm{i})-1,length(rad_lbm{i}));
    eval(sprintf('t_rad_lbm%d = linspace(0,length(rad_lbm%d)-1,length(rad_lbm%d));',i,i,i));
    eval(sprintf('rad_deriv%d = (rad_lbm%d(init_idx+2)-rad_lbm%d(init_idx+1));',i,i,i));
end
cd("matlab_r_p\")
% %%%%%%%%%%% MY-LBM %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% matlab-R-P %%%%%%%%%%%%%%
options = odeset('RelTol', 1e-1, 'AbsTol', 1e-2);
num_files = numel(file_names);
input_files = cell(1, num_files);
% legend_entries3 = cell(num_files, 1);  % Initialize legend entries cell array
for k = 1:num_files
    input_files{1} = sprintf('input%d.json', 27);
    input_files{2} = sprintf('input%d.json', 31);
    paramters2= jsondecode(fileread(input_files{k}));
end

line_styles={'-','--',':','-.','*','square','diamond','x'};
for i = 1:length(input_files)
    eval(sprintf('length_rad_lbm%d = 0;',i));  % Initialize the length to zero
    for l=1:numel(init_idx)
       current_init_idx=init_idx(l)+1;
        for k = current_init_idx:numel(eval(sprintf('rad_lbm%d',i)))
            if eval(sprintf('rad_lbm%d(k)',i)) == 0
                eval(sprintf('length_rad_lbm%d = k - current_init_idx+1;',i));  % Update the length
                break;  % Exit the loop once the condition is met
            else
                eval(sprintf('length_rad_lbm%d = k - current_init_idx+1;',i));
            end
        end
        eval(sprintf('[t, R] = ode45(@(t, R)RP_eval_fct(t,R, input_files{%d}), linspace(t_rad_lbm%d(current_init_idx), length_rad_lbm%d+t_rad_lbm%d(current_init_idx)-1, length_rad_lbm%d), [rad_lbm%d(current_init_idx), rad_deriv%d(l)], options);', i, i, i, i,i,i,i));
        parameters = jsondecode(fileread(input_files{i}));
        r=parameters.radius;
         
        R_{i}{l}=R;
        t_{i}{l}=t;
        hold on;
%         eval(sprintf('difference = rad_lbm%d(current_init_idx:current_init_idx+length_rad_lbm%d-1, 1) - R(1);',i,i));
%         mean_diff=mean(difference);
    legend_entries2{l} = sprintf('R_{R-P,0}=R_{LBM}(t=%d)',init_idx(l));
    legend_entries4{l} = sprintf('|(R_{R-P,0}-R_{LBM}(t=%d)|',init_idx(l));
    end
end
% interpolate
output1 = zeros(1, num_files);  % Preallocate array for output1
output2 = zeros(1, num_files);  % Preallocate array for output2
x=zeros(num_files,numel(init_idx));
y=zeros(num_files,numel(init_idx));
for i= 1:numel(file_names)
    for l=1:numel(init_idx)
%     if(i==myidx)
            if(i==1 && l==1)
                my_lbm_R2{i}{l}=rad_lbm1;
            end
        if ( length(R_{i}{l}(:,1))<length(rad_lbm{i}) )
            my_time = t_{i}{l};
            my_lbm_R = interp1(t_rad_lbm{i},rad_lbm{i},my_time);
            rel_diff = (my_lbm_R-R_{i}{l}(:,1))./R_{i}{l}(:,1);
            my_time2{i}{l}=my_time;
            my_lbm_R2{i}{l}=my_lbm_R;
            
            % Check if the condition is met
            if any(abs(rel_diff * 100) > 5)
                % Calculate the time difference
                output1(i) = my_time(find(abs(rel_diff * 100) > 5, 1)) - my_time(1);
                x(i,l)=my_time(find(abs(rel_diff * 100) > 5, 1));
                y(i,l)=rad_lbm{i}(x(i,l));
                % Calculate the radius
%                 output2(i) = eval(sprintf('abs(rad_lbm%d(my_time(find(abs(rel_diff*100)>5,1)))-rad_lbm%d(1))/rad_lbm%d(1)*100', i,i,i));
            else
                % Condition not met, set output to 1000
                output1(i) = 1000;
                x(i,l)=1000;
                y(i,l)=rad_lbm{i}(x(i,l));
%                 output2(i) = R_{i}(length(R_{i}));
            end
        end
    end
end
screenSize = get(groot,'ScreenSize');
set(gcf,'units','normalized','outerposition',[0 0 1 1])

legend_entries = cell(1, 4); 
current_legend_entries = cell(1, 4);
legend_entries3=cell(1,4);
legend_entries3{1}="";
legend_entries3{2}="";
legend_entries3{3}="";
legend_entries3{4}='5% difference';

% for j = 1:8
%     if mod(j, 2) == 1
%         k = (j + 1) / 2;  % For j=1,3,5,7 growth
%         subplot_index = k;
%         i_start = (j - 1) * 4 + 1;
%         i_end = j * 4;
%         hold off;
%         figure(k)
% Create a new pair of axes inside the current subplot

mylinestyle={'-o','-square','-diamond','->'};
for i=1:num_files
    hold off;
    subplot(2, 2, i);
    handles(i) = subplot(2,2,i);
%     ax2 = axes('position', [0.35 0.6 0.1 0.15]);
%     set(ax2, 'box', 'on');
%     ax2.XTickLabel = [];
%     ax2.YTickLabel = [];
    eval(sprintf('plot(t_rad_lbm%d(1:40:end), rad_lbm%d(1:40:end),"-","MarkerSize",3);', i, i));
%     hold on;
    % Create a local zoom-in figure within each subplot
    hold on;
% Plot data
    for l=1:numel(init_idx)
        if(i==2 && l==4)
            t_{i}{l}=t_{i}{l}(1:310);
            R_{i}{l}=R_{i}{l}(1:310)';
        end
%     if any(i >= i_start & i <= i_end)
        % Compute the corresponding position within the current subplot
        eval(sprintf('plot(t_{%d}{%d}(1:80:end), R_{%d}{%d}(1:80:end,1),mylinestyle{l},"LineWidth",1,"MarkerSize",3);', i,l, i,l));
%         scatter(x(i),y(i), 'filled', 'k','Marker', 'o','SizeData',100);
        % Add labels (a), (b), (c), (d) to the top left corner of each subplot
        text(0.9, 0.95, sprintf('(%c)', 'a' + i - 1), 'Units', 'normalized', 'FontSize', 12);
        % Additional plot customization can be added here if needed      
        set(gca,'FontWeight','bold','FontSize',18);
        xlabel('t [t.u.]');
        ylabel('R [l.u.]');
        
        if (i==1)
            xlim([0,1000])
            ylim([20,110]);
        else
            xlim([0,550])
            ylim([0,35]);
        end
        for m=1:5
            if m==1
                legend_entries{m} = legend_entries1{i};
            else
                legend_entries{m} =  legend_entries2{m-1};
            end
        end
        if i==1
            hlegend=legend(legend_entries, 'Location', 'northwest');
            set(hlegend, 'FontSize', 18, 'FontWeight', 'bold');
        end
    end
    set(gca,'FontWeight','bold','FontSize',18);
end

%%local zoom in 
% Create a local zoom-in figure within the second subplot
ax2 = axes('position', [0.35 0.6 0.1 0.15]);
set(ax2, 'box', 'on');
% Get handles of lines in the current axes
line_handles = findobj(handles(1),'Type','line');

% Duplicate specific lines onto ax2
copyobj(line_handles, ax2);
ax2.XTickLabel = [];
ax2.YTickLabel = [];

% Customize ax2
xlim(ax2, [0, 110]);
ylim(ax2, [30, 35]);    
    
% Create a second local zoom-in figure within the second subplot
ax3 = axes('position', [0.58 0.6 0.1 0.15]);
set(ax3, 'box', 'on');
% Get handles of lines in the current axes
line_handles = findobj(handles(2),'Type','line');

% Duplicate specific lines onto ax3
copyobj(line_handles, ax3);
ax3.XTickLabel = [];
ax3.YTickLabel = [];

% Customize ax3
xlim(ax3, [0, 110]);
ylim(ax3, [28, 31]);  
%%

for i=1:num_files
    hold off;
    subplot(2, 2, i+2);
%     eval(sprintf('plot(t_rad_lbm%d, rad_lbm%d,"--");', i, i));
%     hold on;
% Plot data
    for l=1:numel(init_idx)
        if(i==2 && l==4)
            t_{i}{l}=t_{i}{l}(1:256);
            R_{i}{l}=R_{i}{l}(1:256);
            my_lbm_R2{i}{l}=my_lbm_R2{i}{l}(1:256);
        end
%     if any(i >= i_start & i <= i_end)
        % Compute the corresponding position within the current subplot
        eval(sprintf('plot(t_{%d}{%d}(1:40:end), abs(R_{%d}{%d}(1:40:end,1) - my_lbm_R2{%d}{%d}(1:40:end)), ''%s'',"LineWidth",1,"MarkerSize",3);', i, l, i, l, i,l,mylinestyle{l}));
        hold on;
%         scatter(x(i),y(i), 'filled', 'k','Marker', 'o','SizeData',100);
        % Add labels (a), (b), (c), (d) to the top left corner of each subplot
        text(0.9, 0.95, sprintf('(%c)', 'c' + i - 1), 'Units', 'normalized', 'FontSize', 12);
        % Additional plot customization can be added here if needed      
        set(gca,'FontWeight','bold','FontSize',18);
        xlabel('t [t.u.]');
        ylabel('\Delta R [l.u.]');
        
%         if (i==1)
%             xlim([0,1000])
%             ylim([20,110]);
%         else
%             xlim([0,550])
%             ylim([0,35]);
%         end
        for m=1:4
%             if m==1
%                 legend_entries{m} = legend_entries1{i};
%             else
                legend_entries{m} =  legend_entries4{m};
%             end
        end
        if i==1
            hlegend=legend(legend_entries, 'Location', 'northwest');
            set(hlegend, 'FontSize', 18, 'FontWeight', 'bold');
        end
    end
end
%     end
% end
% saveas(gcf,'..\figures\bubble-curve.png')
% run("matlab_r_p/cbar2.m")
cd ../src