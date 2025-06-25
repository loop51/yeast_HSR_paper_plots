clc;
clear;
close all;

%%
clc;
clear;
close all;
plot_line_width = 2.5;
% get september 12th cell data at 32C
if ispc
    experiment_folder = "C:\Users\mdsaifi\Box\Measurement\Fall2024\september12th";
    % figure_save_folder = "C:\Users\mdsaifi\Box\Research\Semesters\Fall2024\Yeast_Thermal_effect_paper\yeast HSR paper draft\Results\March2025\Growth_at_32C_multiple_cells";
    figure_save_folder = "k:\Yeast_Thermal_effect_paper\yeast HSR paper draft\Results\March2025_v2\may\Growth_at_32C_multiple_cells";
elseif ismac
    experiment_folder = '/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Measurement/Fall2024/september12th';
    figure_save_folder = '/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Research/Semesters/Fall2024/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/may/Growth_at_32C_multiple_cells';
elseif isunix && ~ismac
    experiment_folder = '/media/mdsaifi/Grace Lab/measurement/fall2024/september12th/';
    figure_save_folder = '/media/mdsaifi/Grace Lab/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/March2025_v2/Growth_at_32C_multiple_cells/may/';
end

if ~isfolder(figure_save_folder)
    mkdir(figure_save_folder);
end
load(strcat(experiment_folder, filesep, "iteration1", filesep, "mag_phase_temporal_data_1.mat"));
% load(strcat(experiment_folder, filesep, "mag_phase_temporal_data_1.mat"));

sept_12th_cell_data = temporal_data;

% get sept20th data cell data at 32C
if ispc
    load('C:\Users\mdsaifi\Box\Measurement\Fall2024\September20th\Yeast SP at 32C\mag_phase_temporal_data_1.mat');
elseif ismac
    load('/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Measurement/Fall2024/September20th/Yeast SP at 32C/mag_phase_temporal_data_1.mat')
elseif isunix && ~ismac
    load('/media/mdsaifi/Grace Lab/measurement/fall2024/september20th/Yeast SP at 32C/mag_phase_temporal_data_1.mat');
end
    
sept_20th_cell_data = temporal_data;

if ispc
    load('C:\Users\mdsaifi\Box\Measurement\Fall2024\November\28th Nov Yeast YPD\at32C\mag_phase_temporal_data_1.mat');
elseif ismac
    load('/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Measurement/Fall2024/November/28th Nov Yeast YPD/at32C/mag_phase_temporal_data_1.mat')
elseif isunix && ~ismac
    load('/media/mdsaifi/Grace Lab/measurement/fall2024/November/28th Nov Yeast YPD/at32C/mag_phase_temporal_data_1.mat');
end
    
nov_28th_cell_data_didnt_grow = temporal_data;

if ispc
    load('C:\Users\mdsaifi\Box\Measurement\Fall2024\November\28th Nov Yeast YPD\another_cell\at32C\mag_phase_temporal_data_1.mat');
elseif ismac
    load('/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Measurement/Fall2024/November/28th Nov Yeast YPD/at32C/mag_phase_temporal_data_1.mat')
elseif isunix && ~ismac
    load('/media/mdsaifi/Grace Lab/measurement/fall2024/November/28th Nov Yeast YPD/another_cell/at32C/mag_phase_temporal_data_1.mat');
end
    
nov28th_four_cells = temporal_data;

%% modify the originial data absolute time to relative time from the beginning of each sample.
sept_12th_cell_data = config_temporal_data(sept_12th_cell_data);
sept_20th_cell_data = config_temporal_data(sept_20th_cell_data);
nov_28th_cell_data_didnt_grow = config_temporal_data(nov_28th_cell_data_didnt_grow);


nov28th_four_cells  = config_temporal_data(nov28th_four_cells);

%% create labels for each cell
temp = datetime(sept_12th_cell_data(1).data(1,5), 'convertfrom', 'datenum');
sept_12th_label = strcat(string(temp.Month),'/', string(temp.Day));


temp = datetime(sept_20th_cell_data(1).data(1,5), 'convertfrom', 'datenum');
sept_20th_label = strcat(string(temp.Month),'/', string(temp.Day));


temp = datetime(nov_28th_cell_data_didnt_grow(1).data(1,5), 'convertfrom', 'datenum');
nov_28th_label = strcat(string(temp.Month),'/', string(temp.Day));
% nov_28th_label = strcat(string(temp.Month),'/', string(temp.Day));


%%
ii=4;
upl = 100;
ll = 0;

    % S11m
    figure();
    param = 7;
    param = [0, 0.4471, 0.7412];
    deltaS_vs_time_plot(sept_12th_cell_data(ii), param);
    param = [0.8510, 0.3255, 0.0980];  % Dark yellow
    deltaS_vs_time_plot(sept_20th_cell_data(ii), param);
    
    param = [0.9294, 0.6941, 0.1255]; % orange
    deltaS_vs_time_plot(nov_28th_cell_data_didnt_grow(ii), param);
    param = [0, 0,0];  % black
    deltaS_vs_time_plot(nov28th_four_cells(ii), param);
    % xlim([ll upl]);
    % legend(sept_12th_label, sept_20th_label, strcat(nov_28th_label, '_1'), strcat(nov_28th_label, '_2'), 'Location', 'northeast');
    legend('cell1', 'cell2', 'cell3', 'cell4', 'Location', 'northwest');
    % legend(sept_20th_label, nov_28th_label);
    % ylabel('\Delta |S_{11|}',  'FontAngle', 'italic');
    % title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))

    % 
    current_figure = gcf;
    current_plot = gca;
    current_figure.Units = 'inches';
    % current_figure.Position = [21.3229 0.4583 3.5 3.5];
    current_figure.OuterPosition = [21.3229 0.4583 7 3.5];
    current_figure.MenuBar = 'none';
    current_plot.FontSize = 11;

    filename = strcat('S11m_vs_S11a_multiple_cell', num2str(temporal_data(ii).frequency), "GHz",'.pdf');
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');



    % 

% end
%%
function [temporal_data] = config_temporal_data(temporal_data)
    data = hold_all_freq(temporal_data);
    % find start time and end time
    for ii=1:length(data)
        all_time(ii) = data(ii,2);
    end
    experiment_start_time = min(all_time);
    % experiment_end_time = max(all_time);
    
    experiment_start_time = datetime(experiment_start_time, 'convertfrom', 'datenum');
    % experiment_end_time = datetime(experiment_end_time, 'convertfrom', 'datenum');


    for ii=1:size(temporal_data,2)
        try
            temporal_data(ii).data(:, 6) = minutes(datetime(temporal_data(ii).data(:, 5), 'ConvertFrom',  'datenum') - experiment_start_time); 
        catch
        end
    end
end

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

function deltaS_vs_time_plot(temporal_data, param)

    temporal_data.data = sortrows(temporal_data.data, 6);
    % x = temporal_data.data(:, 6);
    x = temporal_data.data(:,8);
    y =temporal_data.data(:,7);


    scatter(x,y, 5,param, 'filled');
    grid on;
    hold on;

    % [x,y] = return_polyfit_vals(x,y);

    % plot(x,y,'LineWidth', 1.5);
    grid on;
    hold on;
    % xlabel("t (min)");
    xlabel("\Delta \angleS_{11}",  'FontAngle', 'italic');
    ylabel('\Delta |S_{11}|',  'FontAngle', 'italic');

%     current_figure = gcf;
%     current_plot = gca;
%     current_figure.Units = 'inches';
%     current_figure.Position = [21.3229 0.4583 3.5 3.5];
%     current_figure.MenuBar = 'none';
%     current_plot.FontSize = 10;
% % plot_s11m.LineWidth = 2;
% 

end

function [x_fit,y_fit] = return_polyfit_vals(x, y)
    p = polyfit(x, y, 25);
    x_fit = linspace(min(x), max(x), 1000);
    y_fit = polyval(p, x_fit);

    
end


