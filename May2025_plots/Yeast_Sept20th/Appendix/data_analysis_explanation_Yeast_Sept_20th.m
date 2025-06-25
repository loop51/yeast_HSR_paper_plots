clear;
clc;
close all;

%%
if isunix && ~ismac
    figure_save_folder = '/media/mdsaifi/Grace Lab/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/March2025_v2/september20th/data_analysis_explanation';
elseif ispc
    figure_save_folder = 'K:\Yeast_Thermal_effect_paper\yeast HSR paper draft\Results\March2025_v2\september20th\data_analysis_explanation';
elseif ismac
    % figure_save_folder = '/volumes/Grace Lab/Yeast_Thermal_effect_paper/Yeast_HSR_ADS_simualtions/ADS_Simulation_after_spring_break/figures';
end

if ~isfolder(figure_save_folder)
    mkdir(figure_save_folder);
end

if isunix && ~ismac
    data_path = '/mnt/project/rfcells/November 2024 Measurements/19th Nov PS 5um/PS 5um at 32C/iteration1/3.5GHz/set6/';
    % fig_path = '/mnt/project/rfcells/November 2024 Measurements/19th Nov PS 5um/PS 5um at 32C/iteration1/3.5GHz/set6/filter_adjusted_Processed_Data_Plot.fig';
elseif ispc
    data_path = 'Z:\September2024_Yeast_Measurements\September20th\Yeast SP at 32C\3.5GHz\set6';
elseif ismac
    % figure_save_folder = '/volumes/Grace Lab/Yeast_Thermal_effect_paper/Yeast_HSR_ADS_simualtions/ADS_Simulation_after_spring_break/figures';
end

[s11m,s11a,s21m,s21a,t,f] = get_data(data_path);
% load .mat file that holds the detrended line, and s_filtered lines for
% all variables
load(fullfile(data_path, 'detrend_data.mat'));
y_axis_range_s11m = 0.0025;
y_axis_range_s11a = 0.008;
y_axis_range_s21m = 0.0025;
y_axis_range_s21a = 0.025;
frequency = num2str(f/1e9);


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
current_figure.Position = [13.0104  1.0625 7 3.5];
% current_figure.OuterPosition = [21.3229 0.4583 3.5 3.5];
% current_plot.YRuler.Exponent = -2;
current_figure.MenuBar = 'none';
current_plot.FontSize = 11;

filename = strcat('data_analysis_explanation_s11a_', frequency, 'GHz','.pdf');
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
