clear;
clc;
close all;

%%
if isunix && ~ismac
    figure_save_folder = '/media/mdsaifi/Grace Lab/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/March2025_v2/DI_Water_5um_PS/typical_data/f2';
elseif ismac
    % figure_save_folder = '/volumes/Grace Lab/Yeast_Thermal_effect_paper/Yeast_HSR_ADS_simualtions/ADS_Simulation_after_spring_break/figures';
end

if ~isfolder(figure_save_folder)
    mkdir(figure_save_folder);
end

if isunix && ~ismac
    data_path = '/mnt/project/rfcells/November 2024 Measurements/19th Nov PS 5um/PS 5um at 32C/iteration1/10.5GHz/set6/';
    % fig_path = '/mnt/project/rfcells/November 2024 Measurements/19th Nov PS 5um/PS 5um at 32C/iteration1/3.5GHz/set6/filter_adjusted_Processed_Data_Plot.fig';
elseif ismac
    % figure_save_folder = '/volumes/Grace Lab/Yeast_Thermal_effect_paper/Yeast_HSR_ADS_simualtions/ADS_Simulation_after_spring_break/figures';
end

[s11m,s11a,s21m,s21a,t,f] = get_data(data_path);
% load .mat file that holds the detrended line, and s_filtered lines for
% all variables
load(fullfile(data_path, 'detrend_data.mat'));
y_axis_range_s11m = 0.0025;
y_axis_range_s11a = 0.025;
y_axis_range_s21m = 0.0025;
y_axis_range_s21a = 0.025;
frequency = num2str(f/1e9);
%% s11m figure
current_figure = figure();
plot(t,s11m);
grid on;
hold on;
plot(t, s_filtered.s11m);
plot(t, detrended.s11m);
ylim([min(s11m)-y_axis_range_s11m max(s11m)+y_axis_range_s11m])

ylabel('|S_{11}| (dB)', 'FontAngle', 'italic');
xlabel('t (s)', 'FontAngle', 'italic');

% current_figure = gcf;
current_plot = gca;
current_figure.Units = 'inches';
current_figure.Position = [8 3.0625 3.5 3.5];
% current_figure.OuterPosition = [21.3229 0.4583 3.5 3.5];
% current_plot.YRuler.Exponent = -2;
current_figure.MenuBar = 'none';
current_plot.FontSize = 8;

filename = strcat('s11m_', frequency, 'GHz','.pdf');
exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');

%% s11a figure
current_figure = figure();
plot(t,s11a);
grid on;
hold on;
plot(t, s_filtered.s11a);
plot(t, detrended.s11a);
ylim([min(s11a)-y_axis_range_s11a max(s11a)+y_axis_range_s11a])

ylabel('\angleS_{11}', 'FontAngle', 'italic');
xlabel('t (s)', 'FontAngle', 'italic');

% current_figure = gcf;
current_plot = gca;
current_figure.Units = 'inches';
current_figure.Position = [13.0104  3.0625 3.5 3.5];
% current_figure.OuterPosition = [21.3229 0.4583 3.5 3.5];
% current_plot.YRuler.Exponent = -2;
current_figure.MenuBar = 'none';
current_plot.FontSize = 8;

filename = strcat('s11a_', frequency, 'GHz','.pdf');
exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');

%% s21m figure
current_figure = figure();
plot(t,s21m);
grid on;
hold on;
plot(t, s_filtered.s21m);
plot(t, detrended.s21m);
ylim([min(s21m)-y_axis_range_s21m max(s21m)+y_axis_range_s21m])

ylabel('|S_{21}| (dB)', 'FontAngle', 'italic');
xlabel('t (s)', 'FontAngle', 'italic');

% current_figure = gcf;
current_plot = gca;
current_figure.Units = 'inches';
current_figure.Position = [8 10.0625 3.5 3.5];
% current_figure.OuterPosition = [21.3229 0.4583 3.5 3.5];
% current_plot.YRuler.Exponent = -2;
current_figure.MenuBar = 'none';
current_plot.FontSize = 8;
filename = strcat('s21m_', frequency, 'GHz''.pdf');
exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');

%% s21a figure
current_figure = figure();
plot(t,s21a);
grid on;
hold on;
plot(t, s_filtered.s21a);
plot(t, detrended.s21a);
ylim([min(s21a)-y_axis_range_s21a max(s21a)+y_axis_range_s21a])

ylabel('\angleS_{21}', 'FontAngle', 'italic');
xlabel('t (s)', 'FontAngle', 'italic');

% current_figure = gcf;
current_plot = gca;
current_figure.Units = 'inches';
current_figure.Position = [13.0104  10.0625 3.5 3.5];
% current_figure.OuterPosition = [21.3229 0.4583 3.5 3.5];
% current_plot.YRuler.Exponent = -2;
current_figure.MenuBar = 'none';
current_plot.FontSize = 8;

filename = strcat('s21a_', frequency, 'GHz','.pdf');
exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');

%% plot all in a subplot
font_size = 8;
current_figure = figure();
current_figure.Units = 'inches';
current_figure.Position = [1  15.0625 8.27 3.5/2];
% current_figure.OuterPosition = [1  15.0625 7.5 3.5/2];
current_figure.PaperUnits = 'inches';
current_figure.PaperSize = [8.27 3.5/2];
% current_figure.MenuBar = 'none';
subplot(1,4,1);

plot(t,s11m);
grid on;
hold on;
plot(t, s_filtered.s11m);
plot(t, detrended.s11m);
ylim([min(s11m)-y_axis_range_s11m max(s11m)+y_axis_range_s11m])


ylabel('|S_{11}| (dB)', 'FontAngle', 'italic');
xlabel('t (s)', 'FontAngle', 'italic');
current_plot = gca;
current_plot.FontSize = font_size;


subplot(1,4,2);

plot(t,s11a);
grid on;
hold on;
plot(t, s_filtered.s11a);
plot(t, detrended.s11a);
ylim([min(s11a)-y_axis_range_s11a max(s11a)+y_axis_range_s11a])

ylabel('\angleS_{11}', 'FontAngle', 'italic');
xlabel('t (s)', 'FontAngle', 'italic');

current_plot = gca;
current_plot.FontSize = font_size;

subplot(1,4,3);
plot(t,s21m);
grid on;
hold on;
plot(t, s_filtered.s21m);
plot(t, detrended.s21m);
ylim([min(s21m)-y_axis_range_s21m max(s21m)+y_axis_range_s21m])

ylabel('|S_{21}| (dB)', 'FontAngle', 'italic');
xlabel('t (s)', 'FontAngle', 'italic');

% current_figure = gcf;
current_plot = gca;
current_plot.FontSize = font_size;


subplot(1,4,4);
plot(t,s21a);
grid on;
hold on;
plot(t, s_filtered.s21a);
plot(t, detrended.s21a);
ylim([min(s21a)-y_axis_range_s21a max(s21a)+y_axis_range_s21a])

ylabel('\angleS_{21}', 'FontAngle', 'italic');
xlabel('t (s)', 'FontAngle', 'italic');

current_plot = gca;
current_plot.FontSize = font_size;


filename = strcat('all_s_parms_', frequency, 'GHz','.pdf');
exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');

%%
function [s11m,s11a,s21m,s21a, x_axis_time, CW_Freq] = get_data(data_path)
    % cf = pwd;
    all_files = dir(fullfile(data_path, '**', '*'));  % recursively get all files
    all_files = all_files(~[all_files.isdir]);       % remove directories

    % Filter files using regexp
    match_idx = ~cellfun(@isempty, regexp({all_files.name}, 'S11'));
    temp = all_files(match_idx);
    disp(temp.name);
    data_s11 = csvread(fullfile(temp.folder, temp.name));

    match_idx = ~cellfun(@isempty, regexp({all_files.name}, 'S21'));
    temp = all_files(match_idx);
    disp(temp.name);
    data_s21 = csvread(fullfile(temp.folder, temp.name));

    match_idx = ~cellfun(@isempty, regexp({all_files.name}, '_x_axis_time'));
    temp = all_files(match_idx);
    disp(temp.name);
    x_axis_time = csvread(fullfile(temp.folder, temp.name));

    s11m = 20*log10(abs(data_s11)); 
    s11m = 20*log10(abs(data_s11));
    s21m = 20*log10(abs(data_s21));
    s11a = rad2deg(angle(data_s11));
    s21a = rad2deg(angle(data_s21));


    file_name = temp.name;

    parsed = regexp(file_name, '\_', 'split');

    date_stamp = string(parsed(1));

    SenseType = string(parsed(4));

    CW_Freq = parsed(6);
    CW_Freq = str2num(string(CW_Freq)) *1e9;
end
function [nama] = only_get_main_folder(temp)
    cf = pwd;
    more_temp = temp;
    for ii=1:length(more_temp)
        if string(more_temp(ii).folder) == string(cf)
            temp = more_temp(ii);
        end
    end
    nama = temp;
end
