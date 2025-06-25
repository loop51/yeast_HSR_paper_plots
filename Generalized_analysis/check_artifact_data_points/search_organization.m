clc;
clear;

if ispc
    % experiment_folder = "C:\Users\mdsaifi\Box\Measurement\Fall2024\November\19th Nov PS 5um";
    % figure_save_folder = "C:\Users\mdsaifi\Box\Research\Semesters\Fall2024\Yeast_Thermal_effect_paper\yeast HSR paper draft\Results\March2025\yeast_SP8\deltaS_vs_time";
    % experiment_folder = "K:\measurement\fall2024\september20th\Yeast SP ramp 32C to 42C";
    experiment_folder = "k:\Measurement\Fall2024\November\19th Nov PS 5um";
elseif ismac
    experiment_folder = '/media/mdsaifi/Grace Lab/measurement/fall2024/september20th';
    % figure_save_folder = '/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Research/Semesters/Fall2024/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/deltaS_vs_time';
elseif isunix && ~ismac
    experiment_folder = '/media/mdsaifi/Saiful_Lab/CHO_Drug_Effect/June_13th_Single_Cell_CHO_PF_Oligomycin/Oligomycin_test/second_try/iteration1';
end

% file_location = "C:\Users\mdsaifi\Box\Measurement\Fall2024\september12th\iteration1";
file_location = experiment_folder;

load(strcat(file_location, filesep, "mag_phase_temporal_data_1.mat"));

experiment_start_time = min(temporal_data(1).data(:,5));
experiment_start_time = datetime(experiment_start_time, 'ConvertFrom', 'datenum');

for ii=1:size(temporal_data,2)
  temporal_data(ii).data = sortrows(temporal_data(ii).data, 5);
  time_stamps = datetime(temporal_data(ii).data(:,5), 'ConvertFrom', 'datenum');
  temporal_data(ii).data(:,5) = minutes(time_stamps - experiment_start_time);

end