clc;
clear;
close all;
% Get screen size
% screenSize = get(groot,'ScreenSize');
init_idx=100;
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
legend_entries1=cell(1,numel(parameters.mx)*numel(parameters.rho_l));
file_names=cell(1,numel(parameters.mx)*numel(parameters.rho_l));

% myidx = 1,2,3,4;
rad_lbm = {};
t_rad_lbm = {};
R_={};
t_={};
for i = 1:numel(parameters.mx)
    for j = 1:numel(parameters.rho_l)
            file_names{k} = sprintf('testcase2\\pressure_%.0f_%d_%.2f_30_VSM.dat', parameters.lx(i), parameters.mx(i), parameters.rho_l(j));
%             legend_entries1{k} = sprintf('LBM-%d-%.2f-%d', parameters.lx(i), parameters.rho_l(j),parameters.radius(l));
            legend_entries1{k} = sprintf('LBM-$\\rho_{\\ell}$=%.2f',parameters.rho_l(j));
            k=k+1;
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
    eval(sprintf('rad_deriv%d = (rad_lbm%d(init_idx+1)-rad_lbm%d(init_idx-1))/2;',i,i,i));
end
cd("matlab_r_p\")
% %%%%%%%%%%% MY-LBM %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% matlab-R-P %%%%%%%%%%%%%%
options = odeset('RelTol', 1e-1, 'AbsTol', 1e-2);
num_files = numel(file_names);
input_files = cell(1, num_files);
l=1;
for k = 3:4:32
    input_files{l} = sprintf('input%d.json', k);
    paramters2= jsondecode(fileread(input_files{l}));
    l=l+1;
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
    legend_entries2{i} = sprintf('RP-$\\rho_{\\ell}$=%.2f', parameters.rho_l);
    R_{i}=R;
    t_{i}=t;
        hold on;
    eval(sprintf('difference = rad_lbm%d(init_idx:init_idx+length_rad_lbm%d-1, 1) - R(1);',i,i));
    mean_diff=mean(difference);
end
% interpolate
output1 = zeros(1, num_files);  % Preallocate array for output1
output2 = zeros(1, num_files);  % Preallocate array for output2
x=zeros(1,num_files);
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
                output2(i) = eval(sprintf('abs(rad_lbm%d(my_time(find(abs(rel_diff*100)>5,1)))-rad_lbm%d(1))/rad_lbm%d(1)*100', i,i,i));
            else
                % Condition not met, set output to 1000
                output1(i) = 1000;
                x(i)=1000;
                y(i)=rad_lbm{i}(x(i));
                output2(i) = R_{i}(length(R_{i}));
            end
        end
%     end
end
screenSize = get(groot,'ScreenSize');
set(gcf,'units','normalized','outerposition',[0 0 1 1])

legend_entries = cell(1, 4); 
current_legend_entries = cell(1, 4);
legend_entries3=cell(1,2);
legend_entries3{1}="";
legend_entries3{2}='5% difference';
mylinestyle_lbm={'-','-.'};
mylinestyle_rp={'>','o'};
%%
%filter the data
t_rad_lbm1=t_rad_lbm1(rad_lbm1<50 & rad_lbm1>20);
t_rad_lbm2=t_rad_lbm2(rad_lbm2<50 & rad_lbm2>20);
rad_lbm1=rad_lbm1(rad_lbm1<50 & rad_lbm1>20);
rad_lbm2=rad_lbm2(rad_lbm2<50 & rad_lbm2>20);

t_{1}=t_{1}(R_{1}(:,1)>20 & R_{1}(:,1)<50);
t_{2}=t_{2}(R_{2}(:,1)>20 & R_{2}(:,1)<50);
R_{1}=R_{1}(R_{1}(:,1)>20 & R_{1}(:,1)<50);
R_{2}=R_{2}(R_{2}(:,1)>20 & R_{2}(:,1)<50);

t_rad_lbm3=t_rad_lbm3(rad_lbm3<100 & rad_lbm3>20);
t_rad_lbm4=t_rad_lbm4(rad_lbm4<100 & rad_lbm4>20);
rad_lbm3=rad_lbm3(rad_lbm3<100 & rad_lbm3>20);
rad_lbm4=rad_lbm4(rad_lbm4<100 & rad_lbm4>20);

t_{3}=t_{3}(R_{3}(:,1)>20 & R_{3}(:,1)<100);
t_{4}=t_{4}(R_{4}(:,1)>20 & R_{4}(:,1)<100);
R_{3}=R_{3}(R_{3}(:,1)>20 & R_{3}(:,1)<100);
R_{4}=R_{4}(R_{4}(:,1)>20 & R_{4}(:,1)<100);
%%
textAdded = zeros(1, 4); % Assuming there are 4 subplots as per your code
for j = 1:8
    subplot_index = ceil(j / 2);
    % Compute the corresponding position within the current subplot
    subplot(2, 2, subplot_index);
    if mod(j, 2) == 1
        line_style_lbm = mylinestyle_lbm{1}; % Use the first linestyle for odd j
        line_style_rp = mylinestyle_rp{1};
        linewidth=1;
    else
        line_style_lbm = mylinestyle_lbm{2}; % Use the second linestyle for even j
        line_style_rp = mylinestyle_rp{2};
        linewidth=2;
    end
    line_color = eval(sprintf('plot(t_rad_lbm%d(1:30:end), rad_lbm%d(1:30:end),line_style_lbm,"MarkerSize",3,"LineWidth",linewidth);', j, j));
    hold on;
    if subplot_index==1
    eval(sprintf('plot(t_{%d}(1:5:end), R_{%d}(1:5:end,1),line_style_rp,"MarkerSize",3);', j, j));
    else
    eval(sprintf('plot(t_{%d}(1:30:end), R_{%d}(1:30:end,1),line_style_rp,"MarkerSize",3);', j, j));
    end
    scatter(x(j),y(j), 'filled', 'Marker', 'o','SizeData',100, 'MarkerFaceColor', line_color.Color);
% 
%     % Add labels (a), (b), (c), (d) to the top left corner of each subplot
    % Add labels (a), (b), (c), (d) only once per subplot
    if ~textAdded(subplot_index)
        text(0.05, 0.9, sprintf('(%c)', 'a' + subplot_index - 1), 'Units', 'normalized', 'FontSize', 18);
        textAdded(subplot_index) = 1; % Mark as added
    end
    % Additional plot customization can be added here if needed      
    set(gca,'FontSize',18);
    xlabel('t [t.u.]');
    ylabel('R [l.u.]');
    xlim([0,total_time])
    if(subplot_index==1)
        ylim([20,50]);
        xlim([0,1000]);
        yticks([25,30,35,40,45,50]);
    elseif(subplot_index==2)
        ylim([20,100]);
        xlim([0,1000]);
        yticks([40,60,80,100]);
    elseif(subplot_index==3)
        yticks([50,100,150]);
    else
        yticks([20,40,60,80,100,120]);
    end
    
    for l=1:2
        current_legend_entries{1+3*(l-1)} = legend_entries1{l};
        current_legend_entries{2+3*(l-1)} = legend_entries2{l};
        current_legend_entries{3+3*(l-1)} = legend_entries3{l};
    end
    legend_entries= current_legend_entries;
    if subplot_index==3
        hlegend=legend(legend_entries, 'Location', 'northwest');
        set(hlegend, 'FontSize', 14,'Interpreter', 'latex');
    end
end
cd ../src