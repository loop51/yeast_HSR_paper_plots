clc;
clear;
close all

%%
if ispc
    working_folder = 'C:\Users\mdsaifi\Box\Measurement\spring2025\yeast_temperature_experiment_repeat\boradband';
    figure_save_folder = 'K:\Yeast_Thermal_effect_paper\yeast HSR paper draft\Results\March2025\broadband\ypd_air_comparison';
end
if ~isfolder(figure_save_folder)
    mkdir(figure_save_folder);
end
air_folder_location = strcat(working_folder, filesep, 'air', filesep, 'at32C');

s11_air = csvread(strcat(air_folder_location, filesep, '02-Feb-2025_sensor_data_air 32C_IFBW_is_500_NumPoints_is_20001_S11'));
s21_air = csvread(strcat(air_folder_location, filesep, '02-Feb-2025_sensor_data_air 32C_IFBW_is_500_NumPoints_is_20001_S21'));
f_air = csvread(strcat(air_folder_location, filesep, '02-Feb-2025_sensor_data_air 32C_IFBW_is_500_NumPoints_is_20001_Frequency'));

ypd_folder_location = strcat(working_folder, filesep, 'ypd', filesep, 'at32C', filesep, 'set2');

s11_ypd = csvread(strcat(ypd_folder_location, filesep, '04-Feb-2025_sensor_data_YPD 32C_IFBW_is_500_NumPoints_is_20001_S11'));
s21_ypd = csvread(strcat(ypd_folder_location, filesep, '04-Feb-2025_sensor_data_YPD 32C_IFBW_is_500_NumPoints_is_20001_S21'));
f_ypd = csvread(strcat(ypd_folder_location, filesep, '04-Feb-2025_sensor_data_YPD 32C_IFBW_is_500_NumPoints_is_20001_Frequency'));


%%
line_width = 2;
[s11m,s11a] = get_mag_phase(s11_air);
[s21m,s21a] = get_mag_phase(s21_air);
f_air = f_air/1e9;
f_ypd = f_ypd/1e9;
subplot(1,2,1);
hold on;
grid on;
plot(f_air, s11m, 'r-', 'LineWidth', line_width);
current_plot = gca;
current_plot.FontSize = 11;

subplot(1,2,2);
plot(f_ypd, s21m, 'r-','LineWidth', line_width);
current_plot = gca;
current_plot.FontSize = 11;
title('(b)');

subplot(1,2,1);
[s11m,s11a] = get_mag_phase(s11_ypd);
[s21m,s21a] = get_mag_phase(s21_ypd);
title('(a)');
hold on;
grid on;
plot(f_air, s11m, 'k--', 'LineWidth', line_width);
current_plot = gca;
current_plot.FontSize = 11;
xlabel('f (GHz)', 'FontAngle', 'italic');
ylabel('|S_{11}| (dB)', 'FontAngle', 'italic');
legend('air', 'ypd','location', 'southwest');
subplot(1,2,2);
% yyaxis right
hold on;
grid on;
plot(f_air, s21m, 'k--','LineWidth', line_width);

xlabel('f (GHz)', 'FontAngle', 'italic');
ylabel('|S_{21}| (dB)', 'FontAngle', 'italic');

legend('air', 'ypd','location', 'southwest');
% 
current_figure = gcf;
current_plot = gca;
current_plot.FontSize = 11;
current_figure.Units = 'inches';
current_figure.OuterPosition = [21.3229 0.4583 7 3.5];
current_figure.MenuBar = 'none';


% 
filename = strcat('broadband_ypd_air','.pdf');
exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');
filename = strcat('broadband_ypd_air','.png');
exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');

% 

%%
Pin = 10e-3;
i=3493;
Pabs_air = Pin * (1 - abs(s11_air).^2 - abs(s21_air).^2);
Pabs_air(i)
Pabs_ypd = Pin * (1 - abs(s11_ypd).^2 - abs(s21_ypd).^2);
Pabs_ypd(i)
x = Pabs_ypd-Pabs_air;
x(i)
Pabs = x(i);
Cp = 4180;
v = (10e-6) * (100e-6) * (1e-2);
rho = 1000;
m = rho*v;
delta_T =  Pabs/(m*Cp)

function [mag, phase] = get_mag_phase(s)

    mag = 20*log10(abs(s));
    phase = rad2deg(angle(s));
end