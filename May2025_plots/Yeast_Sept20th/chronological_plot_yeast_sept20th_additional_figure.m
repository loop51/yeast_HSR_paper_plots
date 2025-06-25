clear;
clc;
%%
clc;
clear;
close all;
plot_line_width = 2.5;
if ispc
    experiment_folder = "C:\Users\mdsaifi\Box\Measurement\Fall2024\September20th";
    % figure_save_folder = "C:\Users\mdsaifi\Box\Research\Semesters\Fall2024\Yeast_Thermal_effect_paper\yeast HSR paper draft\Results\March2025\yeast_Sept20th\Choronology_Plot";
    figure_save_folder = 'K:\Yeast_Thermal_effect_paper\yeast HSR paper draft\Results\March2025\yeast_Sept20th\Choronology_Plot_all';
elseif ismac

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

%%
for ii=1:size(data_42C, 2)
    temporal_data(ii).data = [data_32C(ii).data; data_42C(ii).data];
end


%% create a data variable that holds all frequency data
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
%% fix color and time
for ii=1:length(data)
    all_time(ii) = data(ii,2);
end
experiment_start_time = min(all_time);
experiment_end_time = max(all_time);

% experiment_start_time=datetime(experiment_start_time, 'ConvertFrom',  'datenum');
% experiment_end_time=datetime(experiment_end_time, 'ConvertFrom',  'datenum');


all_time = sort(all_time);

colors = parula(length(all_time));

for ii=1:length(data)
    current_time = data(ii, 2);
    for jj=1:length(all_time)
        if current_time == all_time(jj)
            data(ii, 7) = colors(jj,1);
            data(ii, 8) = colors(jj,2);
            data(ii, 9) = colors(jj,3);
        end
    end
end

chronology_time = datetime(data(:,2), 'ConvertFrom',  'datenum') - experiment_start_time;



%% use cartesian coordinates
for i=1:length(data)
    % [s11r(i),s11a(i)] = pol2cart(data(i, 3) , data(i,4));
    % [s21m(i),s21m(i)] = pol2cart(data(i, 3) , data(i,4));
    s11m(i) =  data(i, 3);
    s11a(i) =  data(i, 4);

    s21m(i) =  data(i, 5);
    s21a(i) =  data(i, 6);

end

experiment_start_time = datetime(experiment_start_time, 'convertfrom', 'datenum');
experiment_end_time = datetime(experiment_end_time, 'convertfrom', 'datenum');


experiment_end_time = experiment_end_time - experiment_start_time;
% experiment_start_time = experiment_start_time - experiment_start_time;
%% plot s11 with temporal color markings
m = min(datetime(data(:,2), 'ConvertFrom',  'datenum'));
hold on;
% symbols = {'pentagram', '^', 'o', 'square', 'diamond', 'x'};

for ii=1:length(temporal_data)
% ii=4;
    figure();

    for jj=1:length(data)
        if temporal_data(ii).frequency == data(jj, 1)
            hold on;
            elapsed_time = minutes(datetime(data(jj, 2), 'convertfrom', 'datenum') - experiment_start_time);
            if elapsed_time < 50
                scatter(s11a(jj), s11m(jj), 20, data(jj, 7:9), 'filled', 'marker', 'o');
            elseif elapsed_time >50 && elapsed_time < 80
                scatter(s11a(jj), s11m(jj), 20, data(jj, 7:9), 'marker', 'o');
            elseif elapsed_time >=80
                scatter(s11a(jj), s11m(jj), 30, data(jj, 7:9), 'filled', 'marker', 's');
            end
            
            % scatter(s11m(jj), s11a(jj), 50, 'r', 'filled', 'marker', symbols(ii));
            x= datetime(data(jj,2), 'ConvertFrom',  'datenum');
            % text(s11m(jj), s11a(jj), {string(minutes(x-m))});
        end
        
    end
    grid on;


    colorbar;
    colormap(parula);
    caxis([minutes(experiment_start_time - experiment_start_time) minutes(experiment_end_time)]);
    % a.Label.String = 'minute';
    % a.Label = "minute";
    % clabel("minute");

    % title(strcat("Frequency:", num2str(temporal_data(ii).frequency), "GHz"));
    xlabel('\angle S11');
    ylabel('|S11|');
    
    current_figure = gcf;
    current_plot = gca;
    current_figure.Units = 'inches';
    current_figure.Position = [21.3229 0.4583 3.5 3.5];
    current_figure.MenuBar = 'none';
    current_plot.FontSize = 10;
% plot_s11m.LineWidth = 2;
    filename = strcat('reflection_', num2str(temporal_data(ii).frequency), "GHz",'.pdf');
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');
    filename = strcat('reflection_', num2str(temporal_data(ii).frequency), "GHz",'.png');
    exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');
    
end

% 
% for ii=1:length(temporal_data)
% % ii=1;
%     figure();
%     % current_figure = figure();
%     for jj=1:length(data)
%         if temporal_data(ii).frequency == data(jj, 1)
%             scatter(s21a(jj), s21m(jj), 50, data(jj, 7:9), 'filled', 'marker', 'square');
%             hold on;
%             x= datetime(data(jj,2), 'ConvertFrom',  'datenum');
%             % text(s21m(jj), s21m(jj), {string(minutes(x-m))});
%             % scatter(s11m(jj), s11a(jj), 50, 'r', 'filled', 'marker', symbols(ii));
%         end
% 
%     end
%     grid on;
% 
%     colorbar;
%     colormap(parula);
%     caxis([minutes(experiment_start_time - experiment_start_time) minutes(experiment_end_time)]);
%     % a.Label.qString = 'minute';
%     title(strcat("Frequency:", num2str(temporal_data(ii).frequency), "GHz"));
%     xlabel('\angle S21');
%     ylabel('|S21|');
% 
%     current_figure = gcf;
%     current_plot = gca;
%     current_figure.Units = 'inches';
%     current_figure.Position = [21.3229 0.4583 3.5 3.5];
%     current_figure.MenuBar = 'none';
%     current_plot.FontSize = 10;
% % plot_s11m.LineWidth = 2;
%     filename = strcat('transmission_', num2str(temporal_data(ii).frequency), "GHz",'.pdf');
%     exportgraphics(current_figure,strcat(figure_save_folder, filesep, filename), 'ContentType','vector');
% 
% end
% 
% 

