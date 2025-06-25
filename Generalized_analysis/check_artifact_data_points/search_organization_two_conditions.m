clc;
clear;
close all;
plot_line_width = 2.5;
if ispc
    % experiment_folder = "K:\measurement\fall2024\september20th";
    % experiment_folder = "K:\measurement\spring2025\yeast_temperature_experiment_repeat\yeast_SP8";
    % figure_save_folder = "C:\Users\mdsaifi\Box\Research\Semesters\Fall2024\Yeast_Thermal_effect_paper\yeast HSR paper draft\Results\March2025\yeast_SP8\deltaS_vs_time";
    % experiment_folder = "k:\Measurement\Fall2024\November\19th Nov PS 5um";
    experiment_folder = "k:\Measurement\Fall2024\september19th";
elseif ismac
    experiment_folder = '/media/mdsaifi/Grace Lab/measurement/fall2024/september20th';
    % figure_save_folder = '/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Research/Semesters/Fall2024/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/deltaS_vs_time';
elseif isunix && ~ismac
    experiment_folder = '/media/mdsaifi/Grace Lab/measurement/fall2024/september20th';
end
% if ~isfolder(figure_save_folder)
%     mkdir(figure_save_folder);
% end
data_32C_folder = "PS_5um_at_32C";
data_42C_folder = "PS_5um_at_32C_to_42C";
% data_32C_folder = "at32C";
% data_42C_folder = "at42C";
load(strcat(experiment_folder, filesep, data_32C_folder, filesep, "mag_phase_temporal_data_1.mat"));
data_32C = temporal_data;

% clear temporal_data;
load(strcat(experiment_folder, filesep, data_42C_folder, filesep, "mag_phase_temporal_data_1.mat"));
data_42C = temporal_data;
% if ispc
%     load("C:\Users\mdsaifi\Box\Measurement\Fall2024\September20th\two_conditions_mag_phase.mat");
% elseif ismac
%     load('/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Measurement/Fall2024/September20th/two_conditions_mag_phase.mat');
% end

%% create a data variable that holds all frequency data
for ii=1:size(data_42C, 2)
    temporal_data(ii).data = [ data_32C(ii).data; temporal_data(ii).data];
end


experiment_start_time = min(temporal_data(1).data(:,5));
experiment_start_time = datetime(experiment_start_time, 'ConvertFrom', 'datenum');
% ii=4;
ii = 1;
for ii=1:size(temporal_data,2)
  temporal_data(ii).data = sortrows(temporal_data(ii).data, 5);
  time_stamps = datetime(temporal_data(ii).data(:,5), 'ConvertFrom', 'datenum');
  temporal_data(ii).data(:,5) = minutes(time_stamps - experiment_start_time);

end