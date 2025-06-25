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
    experiment_folder = "K:\measurement\fall2024\september20th";
    figure_save_folder = "k:\Yeast_Thermal_effect_paper\yeast HSR paper draft\Results\March2025_v2\september20th\Appendix_deltaS_vs_frequency";
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
assign_iteration_32C = [1, 1];
assign_set_32C = [1, 14];

assign_iteration_42C = [1, 2,2,3];
assign_set_42C = [24, 2,7,10];

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
%%
t=[];
for i=1:size(time_stamps,1)
    t(i) = min(time_stamps(i,:));
end
%%
experiment_start_time = min(data_32C(1).data(:,5));

experiment_end_time = max(data_42C(3).data(:,5));

experiment_start_time = datetime(experiment_start_time, 'convertfrom', 'datenum');
experiment_end_time = datetime(experiment_end_time, 'convertfrom', 'datenum');

% experiment_end_time = experiment_end_time -experiment_start_time;
% experiment_start_time = experiment_start_time - experiment_start_time;

colors = (parula(length(t)));
temp = [];
minutes_for_label = [];
for i=1:size(time_stamps,1)
    temp = datetime(t(i), 'convertfrom', 'datenum');
    minutes_for_label(i)  = minutes(temp - experiment_start_time);

    label_list(i) = string(round(minutes_for_label(i), 1));
end


%%

figure_s11m = figure();
hold on;

plot(f, deltaS.s11m(1,:), '-o',  'LineWidth', plot_line_width);
plot(f, deltaS.s11m(2,:), '-o',  'LineWidth', plot_line_width);

plot(f, deltaS.s11m(3,:), '--o', 'LineWidth', plot_line_width);
plot(f, deltaS.s11m(4,:), '--o', 'LineWidth', plot_line_width);
plot(f, deltaS.s11m(5,:), '--o', 'LineWidth', plot_line_width);
plot(f, deltaS.s11m(6,:), '--o', 'LineWidth', plot_line_width);
grid on;

xlabel('f (GHz)', 'FontAngle', 'italic');
ylabel('\Delta |S_{11}|', 'FontAngle', 'italic');
plot_s11m = gca;
legend(label_list, 'Location','southeast');


figure_s11a = figure();
hold on;
plot(f, deltaS.s11a(1,:), '-o',  'LineWidth', plot_line_width);
plot(f, deltaS.s11a(2,:), '-o',  'LineWidth', plot_line_width);

plot(f, deltaS.s11a(3,:), '--o', 'LineWidth', plot_line_width);
plot(f, deltaS.s11a(4,:), '--o', 'LineWidth', plot_line_width);
plot(f, deltaS.s11a(5,:), '--o', 'LineWidth', plot_line_width);
plot(f, deltaS.s11a(6,:), '--o', 'LineWidth', plot_line_width);
grid on;

xlabel('f (GHz)', 'FontAngle', 'italic');
ylabel('\Delta \angle S_{11}', 'FontAngle', 'italic');
plot_s11a = gca;
legend(label_list, 'Location','northeast');

figure_s21m = figure();
hold on;

plot(f, deltaS.s21m(1,:), '-o',  'LineWidth', plot_line_width);
plot(f, deltaS.s21m(2,:), '-o',  'LineWidth', plot_line_width);

plot(f, deltaS.s21m(3,:), '--o', 'LineWidth', plot_line_width);
plot(f, deltaS.s21m(4,:), '--o', 'LineWidth', plot_line_width);
plot(f, deltaS.s21m(5,:), '--o', 'LineWidth', plot_line_width);
plot(f, deltaS.s21m(6,:), '--o', 'LineWidth', plot_line_width);
grid on;

xlabel('f (GHz)', 'FontAngle', 'italic');
ylabel('\Delta |S_{21}|', 'FontAngle', 'italic');
plot_s21m = gca;
legend(label_list, 'Location','northwest');


figure_s21a = figure();
hold on;
plot(f, deltaS.s21a(1,:), '-o',  'LineWidth', plot_line_width);
plot(f, deltaS.s21a(2,:), '-o',  'LineWidth', plot_line_width);

plot(f, deltaS.s21a(3,:), '--o', 'LineWidth', plot_line_width);
plot(f, deltaS.s21a(4,:), '--o', 'LineWidth', plot_line_width);
plot(f, deltaS.s21a(5,:), '--o', 'LineWidth', plot_line_width);
plot(f, deltaS.s21a(6,:), '--o', 'LineWidth', plot_line_width);
grid on;

xlabel('f (GHz)', 'FontAngle', 'italic');
ylabel('\Delta \angle S_{21}', 'FontAngle', 'italic');
plot_s21a = gca;
legend(label_list, 'Location','northwest');

%%
figure_s11m.Units = 'inches';
figure_s11m.Position = [21.3229 5 3.5 3.5];
figure_s11m.MenuBar = 'none';
plot_s11m.FontSize = 10;
% plot_s11m.LineWidth = 2;

exportgraphics(figure_s11m,strcat(figure_save_folder, filesep, 'figure_s11m.pdf'), 'ContentType','vector');
% 
figure_s11a.Units = 'inches';
figure_s11a.Position = [26 5 3.5 3.5];
figure_s11a.MenuBar = 'none';
plot_s11a.FontSize = 10;
% plot_s11a.LineWidth = 2;
exportgraphics(figure_s11a,strcat(figure_save_folder, filesep, 'figure_s11a.pdf'), 'ContentType','vector');


figure_s21m.Units = 'inches';
figure_s21m.Position = [21.3229 0.4583 3.5 3.5];
figure_s21m.MenuBar = 'none';
plot_s21m.FontSize = 10;
% plot_s21m.LineWidth = 2;
exportgraphics(figure_s21m,strcat(figure_save_folder, filesep, 'figure_s21m.pdf'), 'ContentType','vector');

figure_s21a.Units = 'inches';
figure_s21a.Position = [26 0.4583 3.5 3.5];
figure_s21a.MenuBar = 'none';
plot_s21a.FontSize = 10;
% plot_s21a.LineWidth = 2;
exportgraphics(figure_s21a,strcat(figure_save_folder, filesep, 'figure_s21a.pdf'), 'ContentType','vector');




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

        t(i) = min(temp.t);
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

