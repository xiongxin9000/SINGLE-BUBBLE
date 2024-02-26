clc;
clear;
close all;

% Get screen size
screenSize = get(groot,'ScreenSize');

% Load parameters from JSON file
parameters = jsondecode(fileread('input_matlab.json'));
data_re = load("reference.dat");
T_reference = data_re(:, 2);
rho_reference = data_re(:, 1);
rho_l = parameters.rho_l;

% Load data for rho_liquid and rho_gas
rho_liquid = [];
rho_gas = [];
for i = 1:length(rho_l)
    data = load("testcase2/" + rho_l(i) + ".dat");
    rho_liquid = [rho_liquid, data(end, end)];
    rho_gas = [rho_gas, data(end, 1)];
end

% Combine rho_liquid and rho_gas into one single data array rho
rho_liquid=rho_liquid';
rho_gas=flip(rho_gas)';
rho = [rho_liquid; rho_gas];

T_re = parameters.T;
T = [T_re; T_re(end:-1:1)];

% Plot the data
figure;
plot(rho, T, 'Marker', 'o');
hold on;
plot(rho_reference, T_reference, 'Marker', 'x')
set(gca, 'FontSize', 18);

% Add legend for rho and rho_reference
legend('LBM', 'Reference', 'Location', 'northwest');
xlabel('$\rho\, \bf{[l.u.]}$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('T [l.u.]');
% title('Maxwell Area Construction');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);

saveas(gcf, 'figures2/maxwell.png');
cd src\


