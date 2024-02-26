clc;
clear;
close all;
% Get screen size
% screenSize = get(groot,'ScreenSize');
% init_idx=10;
init_idx=100;
% for i=1:numel(init_idx)
%     disp(init_idx(i))
% end
% % Set figure position
% figure('Position',[1 1 screenSize(3) screenSize(4)]);
% %%%%%%%%%%% MY-LBM %%%%%%%%%%%%
cd("..\")
% Define file names
% file_names = ["testcase2\pressure_100.00_101_0.31_VSM.dat", "testcase1\pressure_100.00_101_0.34_VSM.dat", "testcase1\pressure_400.00_401_0.31_VSM.dat", "testcase1\pressure_400.00_401_0.34_VSM.dat"];
parameters = jsondecode(fileread('input_matlab2.json'));
radius=parameters.radius;
% radius(1)
total_time=parameters.mstep;
% file_names = cell(1, numel(parameters.mx)*numel(parameters.rho_l));
k=1;
file_names=cell(1,5);
legend_entries1=cell(1,numel(file_names));
% myidx = 1,2,3,4;
rad_lbm = {};
t_rad_lbm = {};
R_={};
t_={};
parameters.mx=1000;
parameters.lx=1001;
parameters.rho_l=0.34;
for i = 1:numel(parameters.mx)
    for j = 1:numel(parameters.rho_l)
        for l = 2:numel(parameters.radius)
            file_names{k} = sprintf('testcase2\\pressure_%.0f_%d_%.2f_%d_VSM.dat', parameters.mx(i),parameters.lx(i), parameters.rho_l(j),parameters.radius(l));
            legend_entries1{k} = sprintf('');
            k=k+1;
        end
    end
end
legend_entries1{end-1}='LBM';
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
    input_files{k} = sprintf('input%d.json', 28+k);
    paramters2= jsondecode(fileread(input_files{k}));
end

line_styles={'-','--',':','-.','*','square','diamond','x'};
for i = 1:length(input_files)
    eval(sprintf('length_rad_lbm%d = 0;',i));  % Initialize the length to zero
        for k = init_idx:numel(eval(sprintf('rad_lbm%d',i)))
            if eval(sprintf('rad_lbm%d(k)',i)) == 0
                eval(sprintf('length_rad_lbm%d = k - init_idx+1;',i));  % Update the length
                break;  % Exit the loop once the condition is met
            else
                eval(sprintf('length_rad_lbm%d = k - init_idx+1;',i));
            end
        end
        eval(sprintf('[t, R] = ode45(@(t, R)RP_eval_fct(t,R, input_files{%d}), linspace(t_rad_lbm%d(init_idx), length_rad_lbm%d+t_rad_lbm%d(init_idx)-1, length_rad_lbm%d), [rad_lbm%d(init_idx), rad_deriv%d], options);', i, i, i, i,i,i,i));
        parameters = jsondecode(fileread(input_files{i}));
        r=parameters.radius;
         
        R_{i}=R;
        t_{i}=t;
        hold on;
%         eval(sprintf('difference = rad_lbm%d(current_init_idx:current_init_idx+length_rad_lbm%d-1, 1) - R(1);',i,i));
%         mean_diff=mean(difference);
    legend_entries2{i} = sprintf('');
end
legend_entries2{end-1}='R-P';
% interpolate
output1 = zeros(1, num_files);  % Preallocate array for output1
output2 = zeros(1, num_files);  % Preallocate array for output2
x=zeros(1,num_files);
y=zeros(1,num_files);
for i= 1:numel(file_names)
%     if(i==myidx)
        if ( length(R_{i}(:,1))<length(rad_lbm{i}) )
            my_time = t_{i};
            my_lbm_R = interp1(t_rad_lbm{i},rad_lbm{i},my_time);
            rel_diff = (my_lbm_R-R_{i}(:,1))./R_{i}(:,1);
           
            % Check if the condition is met
            if any(abs(rel_diff * 100) > 5)
                % Calculate the time difference
                output1(i) = my_time(find(abs(rel_diff * 100) > 5, 1)) - my_time(1);
                x(i)=my_time(find(abs(rel_diff * 100) > 5, 1));
                y(i)=rad_lbm{i}(x(i));
                % Calculate the radius
%                 output2(i) = eval(sprintf('abs(rad_lbm%d(my_time(find(abs(rel_diff*100)>5,1)))-rad_lbm%d(1))/rad_lbm%d(1)*100', i,i,i));
            else
                % Condition not met, set output to 1000
                output1(i) = 1000;
                x(i)=1000;
                y(i)=rad_lbm{i}(x(i));
%                 output2(i) = R_{i}(length(R_{i}));
            end
        end
end
screenSize = get(groot,'ScreenSize');
set(gcf,'units','normalized','outerposition',[0 0 1 1])

legend_entries = cell(1, 12); 
current_legend_entries = cell(1, 12);
legend_entries3=cell(1,5);
legend_entries3{1}="";
legend_entries3{2}="";
legend_entries3{3}="";
legend_entries3{4}="";
legend_entries3{5}='5% difference';
mylinestyle_lbm={'-','--','-.',':'};
mylinestyle_lbm2={'-','--','-.',':','-o'};
mylinestyle_rp={'-diamond','-x','-.','--'};
mylinestyle_rp2={'-diamond','-x','-.','--','hexagram'};
for i=1:num_files
    if i~= num_files
       subplot(2,1,1)
       k=1;
        line_color =eval(sprintf('plot(t_rad_lbm%d(rad_lbm%d>0), rad_lbm%d(rad_lbm%d>0),"-k");', i, i,i,i));
        hold on;
        if(i==3)
            t2_{i}=t_{i}(1:350);
            R2_{i}=R_{i}(1:350)';
        else
            t2_{i}=t_{i};
            R2_{i}=R_{i};
        end
        eval(sprintf('plot(t2_{%d}(1:20:end), R2_{%d}(1:20:end,1),"--k");', i, i));
        scatter(x(i),y(i),'filled', 'Marker', 'o','SizeData',100, 'MarkerFaceColor', line_color.Color);
        % Additional plot customization can be added here if needed      
        set(gca,'FontSize',18);
        xlabel('t [t.u.]');
        ylabel('R [l.u.]');
        xlim([0,total_time]);
        ylim([0,40]);
        yticks([10 20 30 40]);
        for l=1:4
            current_legend_entries{1+3*(l-1)} = legend_entries1{l};
            current_legend_entries{2+3*(l-1)} = legend_entries2{l};
            current_legend_entries{3+3*(l-1)} = legend_entries3{l};
        end
%         legend_entries= current_legend_entries;
%         hlegend=legend(legend_entries, 'Location', 'southeast');
%         set(hlegend, 'FontSize', 14);
        if i==1
        text(0.05, 0.95, sprintf('(%c)', 'a' + k - 1), 'Units', 'normalized', 'FontSize', 18);
        end

    end
        %     else
%         subplot(1,3,2)
%         k=2;
%         eval(sprintf('plot(t_rad_lbm%d, rad_lbm%d,"--");', i, i));
%         hold on;
%         eval(sprintf('plot(t_{%d}, R_{%d}(:,1));', i, i));
%         scatter(x(i),y(i),'filled', 'k','Marker', 'o','SizeData',100);
%         % Additional plot customization can be added here if needed      
%         set(gca,'FontWeight','bold','FontSize',18);
%         xlabel('t [l.u.]');
%         ylabel('R [l.u.]');
%         xlim([0,total_time]);
%         ylim([0,400]);
%         for l=1:5
%         current_legend_entries{1+3*(l-1)} = legend_entries1{l};
%         current_legend_entries{2+3*(l-1)} = legend_entries2{l};
%         current_legend_entries{3+3*(l-1)} = legend_entries3{l};
%         end
%         for m=1:3
%             legend_entries{m}= current_legend_entries{end-2+m-1};
%         end
%         hlegend=legend(legend_entries, 'Location', 'southeast');
%         set(hlegend, 'FontSize', 14, 'FontWeight', 'bold'); 
%         text(0.05, 0.95, sprintf('(%c)', 'a' + k - 1), 'Units', 'normalized', 'FontSize', 12);
%     end
    subplot(2,1,2)
    k=2;
    if i~=5
        line_color =eval(sprintf('loglog(t_rad_lbm%d(1:20:end), rad_lbm%d(1:20:end),"-k");', i, i));
    else
        line_color =eval(sprintf('loglog(t_rad_lbm%d([1:20:500 501:80:3514]), rad_lbm%d([1:20:500 501:80:3514]), "-k");', i, i));
    end    
    hold on;
    if(i==3)
            t_{i}=t_{i}(1:350);
            R_{i}=R_{i}(1:350)';
    end
    if i~=5
        eval(sprintf('loglog(t_{%d}(1:40:end), R_{%d}(1:40:end,1),"--k");', i, i));
    else
        eval(sprintf('loglog(t_{%d}(1:80:end), R_{%d}(1:80:end,1),"--k");', i, i));
    end
    scatter(x(i),y(i),'filled', 'Marker', 'o','SizeData',100, 'MarkerFaceColor', line_color.Color);
    % Additional plot customization can be added here if needed      
    set(gca,'FontSize',18);
    xlabel('t [t.u.]');
    ylabel('R [l.u.]');
    xlim([0,3700]);
%     ylim([0,40]);
    for l=1:5
        current_legend_entries{1+3*(l-1)} = legend_entries1{l};
        current_legend_entries{2+3*(l-1)} = legend_entries2{l};
        current_legend_entries{3+3*(l-1)} = legend_entries3{l};
    end
    legend_entries= current_legend_entries;
    hlegend=legend(legend_entries, 'Location', 'southeast',"Orientation","horizontal","NumColumns", 2);
    set(hlegend, 'FontSize', 18); 
    if i==1
    text(0.05, 0.95, sprintf('(%c)', 'a' + k - 1), 'Units', 'normalized', 'FontSize', 18);
    end
end
%     end
% end
% saveas(gcf,'..\figures\bubble-curve.png')
% run("matlab_r_p/cbar2.m")
cd ../src