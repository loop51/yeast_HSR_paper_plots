clc;
clear;
close all;

%%
% load data
% if ispc
%     load('C:\Users\mdsaifi\Box\Measurement\Fall2024\September20th\two_conditions_mag_phase');
% elseif ismac
% end
% 

%%
clc;
clear;
close all;
plot_line_width = 2.5;
if ispc
    % experiment_folder = "C:\Users\mdsaifi\Box\Measurement\Fall2024\November\19th Nov PS 5um";
    experiment_folder = "k:\measurement\fall2024\November\19th Nov PS 5um";
    figure_save_folder = "k:\Yeast_Thermal_effect_paper\yeast HSR paper draft\Results/March2025_v2\DI_Water_5um_PS\Appendix_deltaS_vs_frequency";
elseif ismac
    experiment_folder = '/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Measurement/Fall2024/November/19th Nov PS 5um';
    figure_save_folder = '/Users/ruzvaysaiful/Library/CloudStorage/Box-Box/Research/Semesters/Fall2024/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/March2025_v2/DI_Water_5um_PS/Appendix_deltaS_vs_frequency';
elseif isunix && ~ismac
    experiment_folder = '/media/mdsaifi/Grace Lab/measurement/fall2024/November/19th Nov PS 5um';
    figure_save_folder = '/media/mdsaifi/Grace Lab/Yeast_Thermal_effect_paper/yeast HSR paper draft/Results/March2025_v2/DI_Water_5um_PS/Appendix_deltaS_vs_frequency';
end
if ~isfolder(figure_save_folder)
    mkdir(figure_save_folder);
end
data_32C_folder = "PS 5um at 32C";
data_42C_folder = "PS 5um 42C";
load(strcat(experiment_folder, filesep, data_32C_folder, filesep, "mag_phase_temporal_data_1.mat"));
data_32C = temporal_data;
% clear temporal_data;
load(strcat(experiment_folder, filesep, data_42C_folder, filesep, "mag_phase_temporal_data_1.mat"));
data_42C = temporal_data;


%%
assign_iteration_32C = [1, 1, 1];
assign_set_32C = [1, 10, 20];

assign_iteration_42C = [1, 1, 2, 2];
assign_set_42C = [5, 14, 10, 15];

%%
deltaS = [];
time_stamps = [];
for i=1:length(assign_iteration_32C)
    [s,t,f] = get_mean(data_32C, assign_iteration_32C(i), assign_set_32C(i));
    % plot(f,s.s11m);
    deltaS.s11m(i,:) = s.s11m;
    deltaS.s11a(i,:) = s.s11a;

    deltaS.s21m(i,:) = s.s21a;
    deltaS.s21a(i,:) = s.s21a;

    time_stamps(i,:) = t;
end
k=size(deltaS.s11m,1)+1;

for i=1:length(assign_iteration_42C)
    [s,t,f] = get_mean(data_42C, assign_iteration_42C(i), assign_set_42C(i));
    % plot(f,s.s11m);
    deltaS.s11m(k,:) = s.s11m;
    deltaS.s11a(k,:) = s.s11a;

    deltaS.s21m(k,:) = s.s21a;
    deltaS.s21a(k,:) = s.s21a;

    time_stamps(k,:) = t;
    k=k+1;
end
t=[];
for i=1:size(time_stamps,1)
    t(i) = mean(time_stamps(i,:));

end

experiment_start_time = min(t);
experiment_end_time = max(t);

experiment_start_time = datetime(experiment_start_time, 'convertfrom', 'datenum');
experiment_end_time = datetime(experiment_end_time, 'convertfrom', 'datenum');

experiment_end_time = experiment_end_time -experiment_start_time;
experiment_start_time = experiment_start_time - experiment_start_time;

colors = (parula(length(t)));


%%

figure_s11m = figure();

for i=1:size(deltaS.s11m,1)
    hold on;
    if i<4
        plot(f, deltaS.s11m(i,:), '-o', 'Color', colors(i,:), 'LineWidth', plot_line_width);
    else
        plot(f, deltaS.s11m(i,:), '--o','Color', colors(i,:), 'LineWidth', plot_line_width);
    end
end
grid on;
colorbar;
colormap(parula);
caxis([minutes(experiment_start_time) minutes(experiment_end_time)]);

% ylim([-1e-3 4.5e-3]);
% for i=1:length(t)
    % temp = minutes(datetime(min(t), 'ConvertFrom', 'datenum') - (datetime(t(i), 'ConvertFrom', 'datenum')));
    % legend(strcat(temp, ' min'));
% end
xlabel('f (GHz)', 'FontAngle', 'italic');
ylabel('\Delta |S_{11}|', 'FontAngle', 'italic');

plot_s11m = gca;

figure_s11a = figure();
for i=1:size(deltaS.s11m,1)
    hold on;
    if i<4
        plot(f, deltaS.s11a(i,:), '-o', 'Color', colors(i,:), 'LineWidth', plot_line_width);
    else
        plot(f, deltaS.s11a(i,:), '--o','Color', colors(i,:),'LineWidth', plot_line_width);
    end
end
grid on;
colorbar;
colormap(parula);
caxis([minutes(experiment_start_time) minutes(experiment_end_time)]);
% ylim([-0.01 0.02]);
xlabel('f (GHz)', 'FontAngle', 'italic');
ylabel('\Delta \angle S_{11}', 'FontAngle', 'italic');
plot_s11a = gca;

figure_s21m = figure();
for i=1:size(deltaS.s11m,1)
    hold on;
    if i<4
        plot(f, deltaS.s21m(i,:), '-o', 'Color', colors(i,:), 'LineWidth', plot_line_width);
    else
        plot(f, deltaS.s21m(i,:), '--o','Color', colors(i,:), 'LineWidth', plot_line_width);
    end
end
grid on;
colorbar;
colormap(parula);
caxis([minutes(experiment_start_time) minutes(experiment_end_time)]);
% ylim([-5e-3 0.03]);
xlabel('f (GHz)', 'FontAngle', 'italic');
ylabel('\Delta |S_{21}|', 'FontAngle', 'italic');
plot_s21m = gca;

figure_s21a = figure();
for i=1:size(deltaS.s11m,1)
    hold on;
    if i<4
        plot(f, deltaS.s21a(i,:), '-o', 'Color', colors(i,:), 'LineWidth', plot_line_width);
    else
        plot(f, deltaS.s21a(i,:), '--o','Color', colors(i,:), 'LineWidth', plot_line_width);
    end
end
grid on;
colorbar;
colormap(parula);
caxis([minutes(experiment_start_time) minutes(experiment_end_time)]);
% ylim([-5e-3 0.03]);
xlabel('f (GHz)', 'FontAngle', 'italic');
ylabel('\Delta \angle S_{21}', 'FontAngle', 'italic');
plot_s21a = gca;
%%
figure_dimension = 3.5;
figure_s11m.Units = 'inches';
figure_s11m.Position = [1 0.4583 figure_dimension figure_dimension];
figure_s11m.MenuBar = 'none';
plot_s11m.FontSize = 8;
% plot_s11m.LineWidth = 2;

exportgraphics(figure_s11m,strcat(figure_save_folder, filesep, 'figure_s11m.pdf'), 'ContentType','vector');
% exportgraphics(figure_s11m,strcat(figure_save_folder, filesep, 'figure_s11m.eps'), 'ContentType','vector');
% 
figure_s11a.Units = 'inches';

figure_s11a.Position = [5 0.4583 figure_dimension figure_dimension];
figure_s11a.MenuBar = 'none';
plot_s11a.FontSize = 8;
plot_s11a.YRuler.Exponent = -2;
% plot_s11a.LineWidth = 2;
exportgraphics(figure_s11a,strcat(figure_save_folder, filesep, 'figure_s11a.pdf'), 'ContentType','vector');
% exportgraphics(figure_s11a,strcat(figure_save_folder, filesep, 'figure_s11a.eps'), 'ContentType','vector');


figure_s21m.Units = 'inches';
figure_s21m.Position = [10 0.4583 figure_dimension figure_dimension];
figure_s21m.MenuBar = 'none';
plot_s21m.FontSize = 8;
% plot_s21m.LineWidth = 2;
exportgraphics(figure_s21m,strcat(figure_save_folder, filesep, 'figure_s21m.pdf'), 'ContentType','vector');

figure_s21a.Units = 'inches';
figure_s21a.Position = [15 0.4583 figure_dimension figure_dimension];
figure_s21a.MenuBar = 'none';
plot_s21a.FontSize = 8;
% plot_s21a.LineWidth = 2;
exportgraphics(figure_s21a,strcat(figure_save_folder, filesep, 'figure_s21a.pdf'), 'ContentType','vector');

%% plot all in a subplot
font_size = 8;
current_figure = figure();
current_figure.Units = 'inches';
current_figure.Position = [0  0 8.27 1.75];
% current_figure.OuterPosition = [1  15.0625 7.5 3.5/2];
current_figure.PaperUnits = 'inches';
current_figure.PaperSize = [8.27 3.5/2];
%%
subplot(1,4,1);
for i=1:size(deltaS.s11m,1)
    hold on;
    if i<4
        plot(f, deltaS.s11m(i,:), '-o', 'Color', colors(i,:), 'LineWidth', plot_line_width);
    else
        plot(f, deltaS.s11m(i,:), '--o','Color', colors(i,:), 'LineWidth', plot_line_width);
    end
end
grid on;
colorbar;
colormap(parula);
caxis([minutes(experiment_start_time) minutes(experiment_end_time)]);

xlabel('f (GHz)', 'FontAngle', 'italic');
ylabel('\Delta |S_{11}|', 'FontAngle', 'italic');

plot_s11m = gca;

plot_s11m.FontSize = 10;
%%
subplot(1,4,2);

for i=1:size(deltaS.s11m,1)
    hold on;
    if i<4
        plot(f, deltaS.s11a(i,:), '-o', 'Color', colors(i,:), 'LineWidth', plot_line_width);
    else
        plot(f, deltaS.s11a(i,:), '--o','Color', colors(i,:),'LineWidth', plot_line_width);
    end
end
grid on;
colorbar;
colormap(parula);
caxis([minutes(experiment_start_time) minutes(experiment_end_time)]);
% ylim([-0.01 0.02]);
xlabel('f (GHz)', 'FontAngle', 'italic');
ylabel('\Delta \angle S_{11}', 'FontAngle', 'italic');
plot_s11a = gca;
plot_s11a.FontSize = 10;
%%
subplot(1,4,3);
for i=1:size(deltaS.s11m,1)
    hold on;
    if i<4
        plot(f, deltaS.s21m(i,:), '-o', 'Color', colors(i,:), 'LineWidth', plot_line_width);
    else
        plot(f, deltaS.s21m(i,:), '--o','Color', colors(i,:), 'LineWidth', plot_line_width);
    end
end
grid on;
colorbar;
colormap(parula);
caxis([minutes(experiment_start_time) minutes(experiment_end_time)]);
% ylim([-5e-3 0.03]);
xlabel('f (GHz)', 'FontAngle', 'italic');
ylabel('\Delta |S_{21}|', 'FontAngle', 'italic');
plot_s21m = gca;
plot_s21m.FontSize = 10;
%%
subplot(1,4,4);
for i=1:size(deltaS.s11m,1)
    hold on;
    if i<4
        plot(f, deltaS.s21a(i,:), '-o', 'Color', colors(i,:), 'LineWidth', plot_line_width);
    else
        plot(f, deltaS.s21a(i,:), '--o','Color', colors(i,:), 'LineWidth', plot_line_width);
    end
end
grid on;
colorbar;
colormap(parula);
caxis([minutes(experiment_start_time) minutes(experiment_end_time)]);
% ylim([-5e-3 0.03]);
xlabel('f (GHz)', 'FontAngle', 'italic');
ylabel('\Delta \angle S_{21}', 'FontAngle', 'italic');
plot_s21a = gca;
plot_s21a.FontSize = 10;


exportgraphics(current_figure,strcat(figure_save_folder, filesep, 'all_deltaS_vs_frequency_DI_PS.pdf'), 'ContentType','vector');


%%
% functions

function [s, t,f] = get_mean(data, iteration, set)
    for i=1:size(data,2) % cycle through the frequencies
        temp = [];
        k = 1;
        for j=1:size(data(i).data,1)
            if data(i).data(j,1) == iteration && data(i).data(j,2) == set
                temp.s11m(k) = data(i).data(j,7);
                temp.s11a(k) = data(i).data(j,8);

                temp.s21m(k) = data(i).data(j,9);
                temp.s21a(k) = data(i).data(j,10);
                temp.t(k) = data(i).data(j,5);
                k=k+1;
            end
        end
        s.s11m(i) = mean(temp.s11m);
        s.s11a(i) = mean(temp.s11a);
        
        s.s21m(i) = mean(temp.s21m);
        s.s21a(i) = mean(temp.s21a);

        t(i) = mean(temp.t);
        f(i) = data(i).frequency;
    end
    [s,t,f] = sort_frequency_s_params(s,t,f);
end


function [s, t,f] = sort_frequency_s_params(s,t,f)
    % CUSTOMSORT Sorts an array in ascending order using bubble sort
    % Input: arr - array of numbers to be sorted
    % Output: sorted_array - sorted array in ascending order
    
    % Make a copy of input array
    % f = arr;
    n = length(f);
    
    % Bubble sort implementation
    for i = 1:n-1
        for j = 1:n-i
            % Compare adjacent elements
            if f(j) > f(j+1)
                % Swap elements if they are in wrong order
                temp = f(j);
                f(j) = f(j+1);
                f(j+1) = temp;
                
                temp = s.s11m(j);
                s.s11m(j) = s.s11m(j+1);
                s.s11m(j+1) = temp;

                temp = s.s11a(j);
                s.s11a(j) = s.s11a(j+1);
                s.s11a(j+1) = temp;

                temp = s.s21m(j);
                s.s21m(j) = s.s21m(j+1);
                s.s21m(j+1) = temp;

                temp = s.s21a(j);
                s.s21a(j) = s.s21a(j+1);
                s.s21a(j+1) = temp;

                temp=t(j);
                t(j) = t(j+1);
                t(j+1) = temp;
            end
        end
    end
end

