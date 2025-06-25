clc;
clear;
close all;

%%
clc;
clear;
close all;
plot_line_width = 2.5;
if ispc
    experiment_folder = "C:\Users\mdsaifi\Box\Measurement\spring2025\yeast_temperature_experiment_repeat\yeast_SP8\at42C\iteration1";
    figure_save_folder = "C:\Users\mdsaifi\Box\Research\Semesters\Fall2024\Yeast_Thermal_effect_paper\yeast HSR paper draft\Results\February2025\Yeast_SP8\at42C";
elseif ismac
    experiment_folder = "/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Measurement/spring2025/yeast_temperature_experiment_repeat/yeast_SP8/at32C/iteration1";
    figure_save_folder = '/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Research/Semesters/Fall2024/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/February2025/Yeast_SP8';
end
% Check if the folder exists
if ~isfolder(figure_save_folder)
    % If the folder does not exist, create it
    mkdir(figure_save_folder);
    % disp(['Folder created: ', figure_save_folder]);
else
    % disp(['Folder already exists: ', figure_save_folder]);
end

% data_32C_folder = "Yeast SP at 32C";
% data_42C_folder = "Yeast SP ramp 32C to 42C";
load(strcat(experiment_folder, filesep, "mag_phase_temporal_data_1.mat"));
data_32C = temporal_data;
data_32C = hold_all_freq(temporal_data);
temperature_switch_time = find_last_data_time(data_32C);


%% create a data variable that holds all frequency data
data = hold_all_freq(temporal_data);

%% find start time and end time
for ii=1:length(data)
    all_time(ii) = data(ii,2);
end
experiment_start_time = min(all_time);
experiment_end_time = max(all_time);

experiment_start_time = datetime(experiment_start_time, 'convertfrom', 'datenum');
experiment_end_time = datetime(experiment_end_time, 'convertfrom', 'datenum');

%% set time boxes in the plot
% for_32C_x = linspace(0, minutes(temperature_switch_time-experiment_start_time),1000);
for_32C_x = minutes(temperature_switch_time-experiment_start_time);

% for_transition_x = linspace(minutes(temperature_switch_time-experiment_start_time), (minutes(temperature_switch_time-experiment_start_time)+minutes(30)),1000);

for_transition_x = (minutes(temperature_switch_time-experiment_start_time)+(30));

%%
ii=1;
for ii=1:size(temporal_data,2)
    % S11m
    temporal_data(ii).data = sortrows(temporal_data(ii).data, 5);
    x = minutes(datetime(temporal_data(ii).data(:, 5), 'ConvertFrom',  'datenum') - experiment_start_time);
    y=temporal_data(ii).data(:,7);
    
    deltaS_vs_time_plot(x,y, for_32C_x,for_transition_x);
    ylabel('\Delta |S11|');
    title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))
    
    filename = strcat('S11m_', num2str(temporal_data(ii).frequency), "GHz",'.pdf');
    % xlim([150 220]);
    current_figure = gcf;
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');


    % S11a
    x = minutes(datetime(temporal_data(ii).data(:, 5), 'ConvertFrom',  'datenum') - experiment_start_time);
    y=temporal_data(ii).data(:,8);
    
    deltaS_vs_time_plot(x,y, for_32C_x,for_transition_x);
    ylabel('\Delta \angle S11');
    title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))
    
    filename = strcat('S11a_', num2str(temporal_data(ii).frequency), "GHz",'.pdf');
    % xlim([150 220]);
    current_figure = gcf;
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');


    % S21m
    x = minutes(datetime(temporal_data(ii).data(:, 5), 'ConvertFrom',  'datenum') - experiment_start_time);
    y=temporal_data(ii).data(:,9);
    
    deltaS_vs_time_plot(x,y, for_32C_x,for_transition_x);
    ylabel('\Delta |S21|');
    title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))
    
    filename = strcat('S21m_', num2str(temporal_data(ii).frequency), "GHz",'.pdf');
    % xlim([150 220]);
    current_figure = gcf;
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');


    % S21a
    x = minutes(datetime(temporal_data(ii).data(:, 5), 'ConvertFrom',  'datenum') - experiment_start_time);
    y=temporal_data(ii).data(:,10);
    
    deltaS_vs_time_plot(x,y, for_32C_x,for_transition_x);
    ylabel('\Delta \angle S21');
    title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))
    
    filename = strcat('S21a_', num2str(temporal_data(ii).frequency), "GHz",'.pdf');
    % xlim([150 220]);
    current_figure = gcf;
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');
end


    

%%
function [data] = hold_all_freq(temporal_data)
    kk=1;
    for ii=1:length(temporal_data)
        for jj=1:size(temporal_data(ii).data,1)
            absolute_time = temporal_data(ii).data(jj, 5);
            s11m = temporal_data(ii).data(jj, 7);
            s11a = temporal_data(ii).data(jj, 8);
            s21m = temporal_data(ii).data(jj, 9);
            s21a = temporal_data(ii).data(jj, 10);
            data(kk,:) = [temporal_data(ii).frequency absolute_time s11m s11a s21m s21a];
            kk=kk+1;
        end
    end
end

function t = find_last_data_time(data)
    t = max(data(:,2));
    t= datetime(t, 'convertfrom', 'datenum') + seconds(10);
end

function deltaS_vs_time_plot(data_time,data_S,at_32_time, transition_time)
    x = data_time;
    y = data_S;
    for_32C_x = at_32_time;
    for_transition_x = transition_time;
    figure();

    scatter(x, y, 'filled');
    hold on

    % Create a polynomial fit (change 3 to adjust polynomial degree)
    p = polyfit(x, y, 25);
    x_fit = linspace(min(x), max(x), 1000);
    y_fit = polyval(p, x_fit);
    % plot(x_fit, y_fit, 'r-', 'LineWidth', 1.5)
    xlabel("time in minutes");
    % ylabel('\Delta |S11|');
    % title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))

    %%
    x = x_fit(x_fit<=for_32C_x);
    y = y_fit(1:length(x));
    
%     fill([x fliplr(x)], [y zeros(size(y))], 'b', 'FaceAlpha', 0.3);
% %%
%     x = x_fit(x_fit>for_32C_x & x_fit<=for_transition_x);
%     temp = y;
%     y = y_fit(length(temp)+1:length(temp)+length(x));
% 
%     fill([x fliplr(x)], [y zeros(size(y))], 'y', 'FaceAlpha', 0.3);
% %%
%     x = x_fit(x_fit>for_transition_x);
%     % temp = y;
%     y = y_fit(end-length(x)+1:end);
% 
%     fill([x fliplr(x)], [y zeros(size(y))], 'r', 'FaceAlpha', 0.3);
%     % y_32C = y;


    current_figure = gcf;
    current_plot = gca;
    current_figure.Units = 'inches';
    current_figure.Position = [21.3229 0.4583 3.5 3.5];
    current_figure.MenuBar = 'none';
    current_plot.FontSize = 10;
% plot_s11m.LineWidth = 2;
    

end