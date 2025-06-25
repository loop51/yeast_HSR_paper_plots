clc;
clear;
close all;

%%
clc;
clear;
close all;
plot_line_width = 2.5;
if ispc
    experiment_folder = "C:\Users\mdsaifi\Box\Measurement\Fall2024\September20th";
    figure_save_folder = "k:\Yeast_Thermal_effect_paper\yeast HSR paper draft\Results\March2025_v2\september20th\deltaS_vs_time";
elseif ismac
    experiment_folder = '/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Measurement/Fall2024/September20th';
    figure_save_folder = '/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Research/Semesters/Fall2024/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/deltaS_vs_time';
elseif isunix && ~ismac
    experiment_folder = '/media/mdsaifi/Grace Lab/measurement/fall2024/september20th';
    figure_save_folder = '/media/mdsaifi/Grace Lab/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/March2025_v2/september20th/deltaS_vs_time';
end
if ~isfolder(figure_save_folder)
    mkdir(figure_save_folder);
end
data_32C_folder = "Yeast SP at 32C";
data_42C_folder = "Yeast SP ramp 32C to 42C";
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


data_32C = hold_all_freq(data_32C);
temperature_switch_time = find_last_data_time(data_32C);
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
ii=4;
% for ii=1:size(temporal_data,2)
    % S11m
    x = minutes(datetime(temporal_data(ii).data(:, 5), 'ConvertFrom',  'datenum') - experiment_start_time);
    y=temporal_data(ii).data(:,7);
    
    deltaS_vs_time_plot(x,y, for_32C_x,for_transition_x);
    ylabel('\Delta |S11|');
    % title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))
    current_figure = gcf;

    filename = strcat('S11m_', num2str(temporal_data(ii).frequency), "GHz",'.pdf');
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');

    filename = strcat('S11m_', num2str(temporal_data(ii).frequency), "GHz",'.png');
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');


    % S11a
    x = minutes(datetime(temporal_data(ii).data(:, 5), 'ConvertFrom',  'datenum') - experiment_start_time);
    y=temporal_data(ii).data(:,8);
    
    deltaS_vs_time_plot(x,y, for_32C_x,for_transition_x);
    ylabel('\Delta \angle S11');
    % title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))
    
    filename = strcat('S11a_', num2str(temporal_data(ii).frequency), "GHz",'.pdf');
    current_figure = gcf;
    ax = gca;
    ax.YAxis.Exponent = -2;  
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');
    filename = strcat('S11a_', num2str(temporal_data(ii).frequency), "GHz",'.png');
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');


    % S21m
    x = minutes(datetime(temporal_data(ii).data(:, 5), 'ConvertFrom',  'datenum') - experiment_start_time);
    y=temporal_data(ii).data(:,9);
    
    deltaS_vs_time_plot(x,y, for_32C_x,for_transition_x);
    ylabel('\Delta |S21|');
    % title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))
    
    filename = strcat('S21m_', num2str(temporal_data(ii).frequency), "GHz",'.pdf');
    current_figure = gcf;
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');
    filename = strcat('S21m_', num2str(temporal_data(ii).frequency), "GHz",'.png');
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');

    % S21a
    x = minutes(datetime(temporal_data(ii).data(:, 5), 'ConvertFrom',  'datenum') - experiment_start_time);
    y=temporal_data(ii).data(:,10);
    
    deltaS_vs_time_plot(x,y, for_32C_x,for_transition_x);
    ylabel('\Delta \angle S21');
    % title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))
    
    filename = strcat('S21a_', num2str(temporal_data(ii).frequency), "GHz",'.pdf');
    current_figure = gcf;
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');

    filename = strcat('S21a_', num2str(temporal_data(ii).frequency), "GHz",'.png');
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');
% end


    

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

    scatter(x, y, 5, 'filled');
    hold on

    % Create a polynomial fit (change 3 to adjust polynomial degree)
    p = polyfit(x, y, 25);
    x_fit = linspace(min(x), max(x), 1000);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'r-', 'LineWidth', 1.5)
    xlabel("{\it t (mins)}");
    % ylabel('\Delta |S11|');
    % title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))

    %%
%     x = x_fit(x_fit<=for_32C_x);
%     y = y_fit(1:length(x));
% 
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
    time_points = [0 53 75 86 117.5 184 195 201 327];
    
    for i=1:length(time_points)
        xline(time_points, '-.k');
    end

    x =x_fit(x_fit<=for_32C_x);
    xline(max(x), '--r');

    % text(mean(x)-1.8*std(x), mean(y)-1.75*std(y), '@32C');
    x =x_fit(x_fit<=for_transition_x);
    
    
    xline(max(x), '--r');
    x =x_fit(x_fit>=for_transition_x);
    % text(mean(x), mean(y)-1.75*std(y), '@42C');



    current_figure = gcf;
    current_plot = gca;
    current_figure.Units = 'inches';
    % current_figure.Position = [21.3229 0.4583 3.5 3.5];
    current_figure.OuterPosition = [21.3229 0.4583 3.5 3.5];
    current_figure.MenuBar = 'none';
    current_plot.FontSize = 8;
    % current_plot.Box = 'on';
    % current_figure.OuterPosition
% plot_s11m.LineWidth = 2;
    

end