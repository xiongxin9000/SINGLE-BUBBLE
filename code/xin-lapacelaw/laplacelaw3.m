clc;
clear;
close all;
% Get screen size
screenSize = get(groot,'ScreenSize');
% load file
parameters = jsondecode(fileread('input_matlab.json'));
% file_names = cell(1, numel(parameters.mx)*numel(parameters.rho_l));
k=1;
file_names=cell(1,numel(parameters.radius));
for i = 1:numel(parameters.radius)
        file_names{k} = sprintf('testcase2\\pressure_%.0f_%d_%.2f_%d_VSM.dat', parameters.lx, parameters.mx,parameters.rho_l,parameters.radius(i));
        k=k+1;
end
%inverseR pressure
path = cell(1, numel(file_names));
mylinestyle = {'-', '--', ':'};
subplot(3,1,1);
R_values=[20,25,30];

for i = 1:numel(file_names)
    % Load data
    eval(sprintf('data%d = load(file_names{i});',i))
    eval(sprintf('inverseR%d=data%d(:,2);',i,i))
    eval(sprintf('pressure%d=data%d(:,1);',i,i))
%     eval(sprintf('std%d=std(inverseR%d);',i,i));
    eval(sprintf('normalizedInverseR%d = (inverseR%d - mean(inverseR%d)) / std(inverseR%d);',i,i,i,i));
    t=linspace(1,10000,20000);
    t1=1:5000;
    linestyle=mylinestyle{i};
    mylinewidth=[1,2,3];
    linewidth=mylinewidth(i);
    eval(sprintf('plot(t1,normalizedInverseR%d(1:5000),"MarkerSize",5,"linestyle",linestyle,"Linewidth",linewidth);', i));
    hold on;
    
    if i==1
   text(0.9, 0.9, sprintf('(%c)', 'a' + i - 1), 'Units', 'normalized', 'FontSize', 18);
    end 
    set(gca,'FontSize',18);
   xlabel('$t\, \mathrm{[t.u.^{-1}]}$','Interpreter', 'latex');
    if(i==1)
        ylabel('$(\frac{1}{R} - \overline{R}) / R_{std}$','Interpreter', 'latex');
    end    
%     saveas(gcf,path{i},'png')
    eval(sprintf('last_100_steps_inverseR%d = inverseR%d(end-99:end);',i,i));
    eval(sprintf('last_100_steps_pressure%d = pressure%d(end-99:end);',i,i));
    eval(sprintf('average_inverseR(i) = mean(last_100_steps_inverseR%d);',i));
    eval(sprintf('average_pressure(i) = mean(last_100_steps_pressure%d);',i));
%     xlim([0,1]);
%     ylim([0.8 1.4]);
end
for i = 1:numel(file_names)
    legendEntries{i} = sprintf('R_{init}=%d', R_values(i)); % Assuming R_values(i) exists
end
legend(legendEntries, 'Orientation', 'horizontal', 'Location', 'southeast');
coefficients = polyfit(average_inverseR, average_pressure, 1);
slope = coefficients(1);
hold off;
subplot(3,1,2);
plot(average_inverseR,average_pressure,'Marker','o');
set(gca,'FontSize',18);
text(0.9, 0.9, '(b)', 'Units', 'normalized', 'FontSize', 18,'FontWeight', 'normal');
xlabel('$1/R\, \mathrm{[l.u.^{-1}]}$','Interpreter', 'latex');
ylabel('$\Delta p\,\, \mathrm{[m.u.\cdot\,l.u.^{-1}\cdot t.u.^{-2}]}$', 'Interpreter', 'latex', 'FontSize', 20);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
% saveas(gcf,'figures2/laplacelaw-mat.png')
% figure(2);
subplot(3,1,3);
% Load parameters from JSON file
parameters2 = jsondecode(fileread('input_matlab2.json'));
data_re = load("reference.dat");
T_reference = data_re(:, 2);
T_reference=T_reference./0.09433;
rho_reference = data_re(:, 1);
rho_reference=rho_reference./0.13044;
rho_l = parameters2.rho_l;

% Load data for rho_liquid and rho_gas
rho_liquid = [];
rho_gas = [];
for i = 1:length(rho_l)
    data = load("testcase3/" + rho_l(i) + ".dat");
    rho_liquid = [rho_liquid, data(end, end)];
    rho_gas = [rho_gas, data(end, 1)];
end

% Combine rho_liquid and rho_gas into one single data array rho
rho_liquid=rho_liquid';
rho_gas=flip(rho_gas)';
rho = [rho_liquid; rho_gas];
rho=rho./0.13044;

T_re = parameters2.T;
T = [T_re; T_re(end:-1:1)];
T=T./0.09433;

% Plot the data
plot(rho, T);
hold on;
plot(rho_reference, T_reference, 'x','MarkerSize',7)
set(gca, 'FontSize', 18);
% text(0.9, 0.9, '(c)', 'Units', 'normalized', 'FontSize', 18);
% Add legend for rho and rho_reference
legend('LBM', 'Reference', 'Location', 'northeast');
xlabel('$\rho$/$\rho_c$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$T/T_c$','Interpreter', 'latex');
% title('Maxwell Area Construction');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
cd 'src'