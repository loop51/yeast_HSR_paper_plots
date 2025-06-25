clc;
clear;
close all;
%%
clc;
clear;
close all;
plot_line_width = 2.5;
if ispc
    experiment_folder = "C:\Users\mdsaifi\Box\Measurement\Fall2024\November\19th Nov PS 5um";
    figure_save_folder = "C:\Users\mdsaifi\Box\Research\Semesters\Fall2024\Yeast_Thermal_effect_paper\yeast HSR paper draft\Results\March2025\DI_Water_5um_PS\deltaS_vs_time";
elseif ismac
    experiment_folder = '/Volumes/Grace Lab/measurement/fall2024/November/19th Nov PS 5um';
    figure_save_folder = '/Volumes/Grace Lab/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/March2025_v2/DI_YPD_5um_PS/deltaS_vs_time';
elseif isunix && ~ismac
    experiment_folder = '/media/mdsaifi/Grace Lab/measurement/fall2024/November/19th Nov PS 5um';
    figure_save_folder = '/media/mdsaifi/Grace Lab/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/March2025_v2/DI_YPD_5um_PS/deltaS_vs_time';
end
if ~isfolder(figure_save_folder)
    mkdir(figure_save_folder);
end
data_32C_folder = "PS 5um at 32C";
data_42C_folder = 'PS 5um 42C';
% data_42C_folder = "Yeast SP ramp 32C to 42C";
load(strcat(experiment_folder, filesep, data_32C_folder, filesep, "mag_phase_temporal_data_1.mat"));
di_data_32C = temporal_data;


clear temporal_data;
load(strcat(experiment_folder, filesep, data_42C_folder, filesep, "mag_phase_temporal_data_1.mat"));
data_42C = temporal_data;

for ii=1:size(data_42C, 2)
    temporal_data(ii).data = [ di_data_32C(ii).data; temporal_data(ii).data];
end

ps_DI = temporal_data;

[di_experiment_start_time, di_for_32C_x, di_for_transition_x] = ...
    get_transition_times(di_data_32C, ps_DI);

%% get ps_YPD data
if ispc
    % experiment_folder = "C:\Users\mdsaifi\Box\Measurement\Fall2024\November\19th Nov PS 5um";
    % figure_save_folder = "C:\Users\mdsaifi\Box\Research\Semesters\Fall2024\Yeast_Thermal_effect_paper\yeast HSR paper draft\Results\March2025\DI_Water_5um_PS\deltaS_vs_time";
elseif ismac
    experiment_folder = '/Volumes/Grace Lab/measurement/fall2024/september19th';
    % figure_save_folder = '/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Research/Semesters/Fall2024/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/DI_Water_5um_PS/deltaS_vs_time';
elseif isunix && ~ismac
    experiment_folder = '/media/mdsaifi/Grace Lab/measurement/fall2024/september19th';
    % figure_save_folder = '/media/mdsaifi/Grace Lab/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/March2025_v2/DI_Water_5um_PS/deltaS_vs_time';
end

data_32C_folder = "PS_5um_at_32C";
data_42C_folder = 'PS_5um_at_32C_to_42C';
% data_42C_folder = "Yeast SP ramp 32C to 42C";
load(strcat(experiment_folder, filesep, data_32C_folder, filesep, "mag_phase_temporal_data_1.mat"));
ypd_data_32C = temporal_data;


clear temporal_data;
load(strcat(experiment_folder, filesep, data_42C_folder, filesep, "mag_phase_temporal_data_1.mat"));
ypd_data_42C = temporal_data;

for ii=1:size(data_42C, 2)
    temporal_data(ii).data = [ ypd_data_32C(ii).data; temporal_data(ii).data];
end

ps_YPD = temporal_data;


[ypd_experiment_start_time, ypd_for_32C_x, ypd_for_transition_x] = ...
    get_transition_times(ypd_data_32C, ps_YPD);
%%

for ii=1:size(ps_YPD,2)
    ps_YPD(ii).data = sortrows(ps_YPD(ii).data, 5);
    ps_DI(ii).data = sortrows(ps_DI(ii).data, 5);
end
%%
ii=4;
% for ii=1:size(temporal_data,2)
    % S11m
    param = 7; % for s11m
    
    x = minutes(datetime(ps_YPD(ii).data(:, 5), 'ConvertFrom',  'datenum') - ypd_experiment_start_time);
    y=ps_YPD(ii).data(:,param);
   
    scatter(x,y,25, 'k','filled');
    hold on;
    xline(ypd_for_32C_x, '--k');
    xline(ypd_for_transition_x, '--k');

    x = minutes(datetime(ps_DI(ii).data(:, 5), 'ConvertFrom',  'datenum') - di_experiment_start_time);
    y=ps_DI(ii).data(:,param);

    scatter(x,y,25, 'r','filled');

    xline(di_for_32C_x, '-.r');

    xline(di_for_transition_x, '-.r');

    grid on;
    
    ylabel('\Delta |S11|');
    title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))
    ylim([0 8e-3]);
    filename = strcat('S11m_', num2str(temporal_data(ii).frequency), "GHz",'.pdf');

    current_figure = gcf;
    current_plot = gca;
    current_figure.Units = 'inches';
    current_figure.Position = [21.3229 0.4583 3.5 3.5];
    current_figure.MenuBar = 'none';
    current_plot.FontSize = 10;
    current_plot.FontName = 'Times New Roman';
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector','Resolution', 300);
    % exportgraphics(current_figure, strcat(figure_save_folder, filesep, filename, '.svg'), ...
    % 'ContentType', 'vector');

  %%
    param = 8; % for s11a
    figure();
    x = minutes(datetime(ps_YPD(ii).data(:, 5), 'ConvertFrom',  'datenum') - ypd_experiment_start_time);
    y=ps_YPD(ii).data(:,param);
   
    h1 = scatter(x,y,25, 'k','filled');
    hold on;
    xline(ypd_for_32C_x, '--k');
    xline(ypd_for_transition_x, '--k');

    x = minutes(datetime(ps_DI(ii).data(:, 5), 'ConvertFrom',  'datenum') - di_experiment_start_time);
    y=ps_DI(ii).data(:,param);

    h2 = scatter(x,y,25, 'r','filled');

    xline(di_for_32C_x, '-.r');

    xline(di_for_transition_x, '-.r');

    grid on;

    legend([h1, h2], {'PS in YPD', 'PS in DI Water'}, 'Location', 'northeast')

    ylabel('\angle S11');
    title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))
    ylim([0 0.055]);
    filename = strcat('S11a_', num2str(temporal_data(ii).frequency), "GHz",'.pdf');
    
    current_figure = gcf;
    current_plot = gca;
    current_figure.Units = 'inches';
    current_figure.Position = [21.3229 0.4583 3.5 3.5];
    current_figure.MenuBar = 'none';
    current_plot.FontSize = 10;
    
    current_plot.FontName = 'Times New Roman';

    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector', 'Resolution', 300);
    % exportgraphics(current_figure, strcat(figure_save_folder, filesep, filename, '.svg'), ...
    % 'ContentType', 'vector');


    %%
    param = 9; % for s21m
    figure();
    x = minutes(datetime(ps_YPD(ii).data(:, 5), 'ConvertFrom',  'datenum') - ypd_experiment_start_time);
    y=ps_YPD(ii).data(:,param);
   
    scatter(x,y,25, 'k','filled');
    hold on;
    xline(ypd_for_32C_x, '--k');
    xline(ypd_for_transition_x, '--k');

    x = minutes(datetime(ps_DI(ii).data(:, 5), 'ConvertFrom',  'datenum') - di_experiment_start_time);
    y=ps_DI(ii).data(:,param);

    scatter(x,y,25, 'r','filled');

    xline(di_for_32C_x, '-.r');

    xline(di_for_transition_x, '-.r');

    grid on;

    ylabel('\Delta |S21|');
    title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))
    % ylim([0 5e-3]);
    filename = strcat('S21m_', num2str(temporal_data(ii).frequency), "GHz",'.pdf');

    current_figure = gcf;
    current_plot = gca;
    current_figure.Units = 'inches';
    current_figure.Position = [21.3229 0.4583 3.5 3.5];
    current_figure.MenuBar = 'none';
    current_plot.FontSize = 10;
    current_plot.FontName = 'Times New Roman';
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector', 'Resolution', 300);
    % exportgraphics(current_figure, strcat(figure_save_folder, filesep, filename, '.svg'), ...
    % 'ContentType', 'vector');

    %%
    param = 10; % for s21a
    figure();
    x = minutes(datetime(ps_YPD(ii).data(:, 5), 'ConvertFrom',  'datenum') - ypd_experiment_start_time);
    y=ps_YPD(ii).data(:,param);
   
    scatter(x,y,25, 'k','filled');
    hold on;
    xline(ypd_for_32C_x, '--k');
    xline(ypd_for_transition_x, '--k');

    x = minutes(datetime(ps_DI(ii).data(:, 5), 'ConvertFrom',  'datenum') - di_experiment_start_time);
    y=ps_DI(ii).data(:,param);

    scatter(x,y,25, 'r','filled');

    xline(di_for_32C_x, '-.r');

    xline(di_for_transition_x, '-.r');

    grid on;

    ylabel('\angle S21');
    title(strcat("@", num2str(temporal_data(ii).frequency),'GHz'))
    % ylim([0 5e-3]);
    filename = strcat('S21a_', num2str(temporal_data(ii).frequency), "GHz",'.pdf');

    current_figure = gcf;
    current_plot = gca;
    current_figure.Units = 'inches';
    current_figure.Position = [21.3229 0.4583 3.5 3.5];
    current_figure.MenuBar = 'none';
    current_plot.FontSize = 10;
    current_plot.FontName = 'Times New Roman';
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector', 'Resolution', 300);
    % exportgraphics(current_figure, strcat(figure_save_folder, filesep, filename, '.svg'), ...
    % 'ContentType', 'vector');


%%

function [experiment_start_time, for_32C_x,for_transition_x] = get_transition_times(data_32C, all_data)

    data_32C = hold_all_freq(data_32C);
    temperature_switch_time = find_last_data_time(data_32C);

    all_data = hold_all_freq(all_data);
    for ii=1:length(all_data)
        all_time(ii) = all_data(ii,2);
    end

    experiment_start_time = min(all_time);

    experiment_end_time = max(all_time);

    experiment_start_time = datetime(experiment_start_time, 'convertfrom', 'datenum');
    experiment_end_time = datetime(experiment_end_time, 'convertfrom', 'datenum');

    % when the temperature is switched from 32C to 42C
    for_32C_x = minutes(temperature_switch_time-experiment_start_time);

    % estimated time when the temperature reaches 42C
    for_transition_x = (minutes(temperature_switch_time-experiment_start_time)+(30));
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

function deltaS_vs_time_plot(data_time,data_S,at_32_time, transition_time)
    x = data_time;
    y = data_S;
    for_32C_x = at_32_time;
    for_transition_x = transition_time;
    figure();

    scatter(x, y, 'filled');
    hold on

    % Create a polynomial fit (change 3 to adjust polynomial degree)
    p = polyfit(x, y, 25); % previously it was 25


    x_fit = linspace(min(x), max(x), 1000);
    y_fit = polyval(p, x_fit);
    % plot(x_fit, y_fit, 'r-', 'LineWidth', 1.5)
    xlabel("time in minutes");
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

    x =x_fit(x_fit<=for_32C_x);
    xline(max(x), '-.r');

    text(mean(x)-1.8*std(x), mean(y)-5*std(y), '@32C');
    x =x_fit(x_fit<=for_transition_x);
    
    
    xline(max(x), '-.r');
    x =x_fit(x_fit>=for_transition_x);
    text(mean(x), mean(y)-5*std(y), '@42C');
    

    current_figure = gcf;
    current_plot = gca;
    current_figure.Units = 'inches';
    current_figure.Position = [21.3229 0.4583 3.5 3.5];
    current_figure.MenuBar = 'none';
    current_plot.FontSize = 10;
% plot_s11m.LineWidth = 2;
    

end