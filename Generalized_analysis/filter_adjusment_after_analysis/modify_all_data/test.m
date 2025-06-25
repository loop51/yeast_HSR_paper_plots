%% 
% this script goes into the data, modifies the sivotzky golay filter hyper parameters
% then saves a plot called filter_adjusted_Processed_Data_Plot.fig
% creates detrend_data.mat --> this holds the detrended data original
% filtered data 
% it then modifies the processed_data.mat delta S values to adjusted
% filtered data so, a good estimation of the segmented signal can be made

% filter_adjusted_Processed_Data_Plot.fig -->
% 1) original signal
% 2) detrended signal
% 3) peak points
% 4) delta S values
% 5) uses four suplots

% detrend_data.mat
% 1) s_filtered
% 2) detrended
% 3) t

% works from the folder where all the frequencies are listed
% 0.5GHz
    % set1
    % set2
    % set3
    % set4
% 1.5GHz
% 3.5GHz
% 4.5GHz
% 5.7GHz
% 10.5GHz

% goes into each frequency a folder and modifies all the sets
% generates-->
% detrend_data.mat
% filter_adjusted_Processed_Data_Plot.fig


%%
clear;
clc;
close all;

cf = pwd;


%% filter params
order = 10;
framelen = 81;


frequnecy_folders = get_folders_list(cf);

for ii=1:length(frequnecy_folders)
    working_folder = strcat(cf, filesep, frequnecy_folders(ii).name);
    disp(strcat("accessing: ",working_folder));
    filter_adjusment(working_folder, order, framelen);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%functions used in this code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter_adjusment(dl, order, framelen)
% is the main function-- it calls other functions for this work



function [dir_list] = get_folders_list(directory_location)
    sf_listing = dir(directory_location);
    j=1;
    for i=1:length(sf_listing)
        if sf_listing(i).name ~= "." && sf_listing(i).name ~= ".." && sf_listing(i).name ~= ".snapshot"
            if sf_listing(i).isdir
                dir_list(j) = sf_listing(i);
                j= j+1;
            end
        end
    end
end
function [pd] = to_cartesian(processedData, detrended, s, t)
    k = 1;
    pd = [];
    for ii=1:size(processedData,1)
        for jj=1:length(t)
            if t(jj) == processedData(ii,5)
                d_s11m = 10^(detrended.s11m(jj)/20);
                d_s21m = 10^(detrended.s21m(jj)/20);
                [d_s11r d_s11i] = pol2cart(detrended.s21a(jj), d_s11m);
                [d_s21r d_s21i] = pol2cart(detrended.s21a(jj), d_s21m);

                s11m = 10^(s.s11m(jj)/20);
                s21m = 10^(s.s21m(jj)/20);
                [s11r s11i] = pol2cart(s.s21a(jj), s11m);
                [s21r s21i] = pol2cart(s.s21a(jj), s21m);
                

                pd(k,:) = [processedData(ii,1) processedData(ii,2) processedData(ii,3)...
                    processedData(ii,4) processedData(ii,5)...
                    (s11r-d_s11r) (s11i-d_s11i) (s21r-d_s21r) (s21i-d_s21i)...
                    processedData(ii,10) processedData(ii,11)];

                k=k+1;
            end
        end
    end
end
%%
function [data] = compile_all_sets(dl)
    dl_list = get_folders_list(dl);
    data = [];
    for i=1:length(dl_list)
        if ispc
            working_folder = strcat(dl, "\", dl_list(i).name);
            try
            load(strcat(working_folder, "\", "processed_data.mat"));
            load(strcat(working_folder, "\", "detrend_data.mat"));
            catch
                disp(strcat(working_folder, "-->doesnt have processsed data"));
            end
        elseif ismac
            working_folder = strcat(dl, "/", dl_list(i).name);
            try
            load(strcat(working_folder, "/", "processed_data.mat"));
            load(strcat(working_folder, "/", "detrend_data.mat"));
            catch
                disp(strcat(working_folder, "-->doesnt have processsed data"));
            end
        end
        processedData = to_cartesian(processedData, detrended, s, t);
        data = [data; processedData];
    end
end
%%
function [f] = get_frequency(frequency_directory)
    x = frequency_directory;
    if ismac
        x = regexp(x, '/', 'split');
        filename_final_data = string(strcat(x(length(x))));
    elseif ispc
        x = regexp(x, '\', 'split');
        filename_final_data = string(strcat(x(length(x))));
    end
    f = str2double(regexp(filename_final_data, '\d+(\.\d+)?', 'match'));
end
%%
function [i] = get_iteration(iteration_folder)
    x = iteration_folder;
    if ismac
        x = regexp(x, '/', 'split');
        filename_final_data = string(strcat(x(length(x))));
    elseif ispc
        x = regexp(x, '\', 'split');
        filename_final_data = string(strcat(x(length(x))));
    end
    i = str2num(regexp(filename_final_data, '\d+(\.\d+)?', 'match'));
end
%% functions that create detrend data

function [nama] = only_get_main_folder(temp, cf)
    % cf = pwd;
    more_temp = temp;
    for ii=1:length(more_temp)
        if string(more_temp(ii).folder) == string(cf)
            temp = more_temp(ii);
        end
    end
    nama = temp;
end

function [s, detrend, t] = get_detrend_from_plot(dl)
    openfig(strcat(dl, filesep, "Processed_Data_Plot.fig"));

    y = subplot(2,2,1);
    g = findobj(y,'type', 'line');
    t = g(3).XData;
    s.s11m = g(4).YData;
    detrend.s11m = g(3).YData;

    y = subplot(2,2,2);
    g = findobj(y,'type', 'line');
    t=g(3).XData;
    s.s11a = g(4).YData;
    detrend.s11a = g(3).YData;

    y = subplot(2,2,3);
    g = findobj(y,'type', 'line');
    t=g(3).XData;
    s.s21m = g(4).YData;
    detrend.s21m = g(3).YData;

    y = subplot(2,2,4);
    g = findobj(y,'type', 'line');
    t=g(3).XData;
    s.s21a = g(4).YData;
    detrend.s21a = g(3).YData;

    close all;
end


function [s] = get_s_params_from_csv(dl)
    % cf = pwd;
    temp =  dir(strcat(dl, filesep, '**/*S11'));
    temp = only_get_main_folder(temp, dl);
    
    data_s11 = csvread(strcat(temp.folder, filesep, temp.name));
    temp =  dir(strcat(dl, filesep, '**/*S21'));
    temp = only_get_main_folder(temp,dl);
    data_s21 = csvread(strcat(temp.folder, filesep, temp.name));

    s.s11m = 20*log10(abs(data_s11));
    s.s21m = 20*log10(abs(data_s21));
    s.s11a = rad2deg(angle(data_s11));
    s.s21a = rad2deg(angle(data_s21));
end

function [s_filtred] = apply_filter(s, order, framelen)
    s_filtred.s11m = sgolayfilt(s.s11m, order, framelen);
    s_filtred.s21m = sgolayfilt(s.s21m, order, framelen);
    s_filtred.s11a = sgolayfilt(s.s11a, order, framelen);
    s_filtred.s21a = sgolayfilt(s.s21a, order, framelen);
end
%%
function filter_adjusment(dl, order, framelen) % here dl is folder/path where all the sets are
    dl_list = get_folders_list(dl);
    % data = [];
    for i=1:length(dl_list)
        working_folder = strcat(dl, filesep, dl_list(i).name);
        disp(strcat("in set: ", dl_list(i).name));
        try
            [s_filtered, detrended,t] = get_detrend_from_plot(working_folder); % fetch s_filtered  detrend_data and t from the processed_data plot
            
            load(strcat(working_folder, filesep, 'processed_data.mat'));
        
        
    
            [s] = get_s_params_from_csv(working_folder);
            [s_filtered] = apply_filter(s, order, framelen);

            save(strcat(working_folder, filesep, "detrend_data",'.mat'),'s_filtered', 'detrended','t'); % save these variables as .mat
            
            T = processedData(:, 5);
            peaks = [];
            k = 1;
            % recompute the delta values
            for ii=1:length(t)
                for jj=1:length(T)
                    if T(jj) == t(ii)
                        peaks.s11m(k) = s_filtered.s11m(ii);
                        peaks.s11a(k) = s_filtered.s11a(ii);
                        peaks.s21m(k) = s_filtered.s21m(ii);
                        peaks.s21a(k) = s_filtered.s21a(ii);
                
                        dS.s11m(k) = s_filtered.s11m(ii) - detrended.s11m(ii);
                        dS.s11a(k) = s_filtered.s11a(ii) - detrended.s11a(ii);
                        dS.s21m(k) = s_filtered.s21m(ii) - detrended.s21m(ii);
                        dS.s21a(k) = s_filtered.s21a(ii) - detrended.s21a(ii);
                        k=k+1;
                    end
                end
            end
            % here update the processData 
            for ii=1:size(processedData,1)
                processedData(ii, 6) = dS.s11m(ii);
                processedData(ii, 7) = dS.s11a(ii); 
                processedData(ii, 8) = dS.s21m(ii);
                processedData(ii, 9) = dS.s21a(ii);
            end
    
            save(strcat(working_folder, filesep, "processed_data",'.mat'),'processedData');
    
            x_axis_time = t;
            %% plot s11m with annotations
            subplot(2,2,1);
            % plot(x_axis_time, s11m,'o');
            plot(x_axis_time, s.s11m);
            hold on;
            plot(x_axis_time, s_filtered.s11m);
            hold on;
            for i = 1:length(dS.s11m)
                text((T(i)), peaks.s11m(i), {num2str(dS.s11m(i)), num2str((T(i)))});
            end
            
            grid on;
            grid minor;
            % plot(x_axis_time, baseline.s11m, 'y');
            plot(x_axis_time, detrended.s11m);
            % plot(T, peaks.s11m, 'o');
            
            % plot(valleySL.t, valleySL.s11m, 'ko');
            % plot(valleySR.t, valleySR.s11m, 'ko');
            
            % plot(x_axis_time, bn.s11m, ' r');
            
            ylabel('|S11| in dB');
            xlabel('time in seconds');
            
            %% plot s11a with annotations
            subplot(2,2,2);
            % plot(x_axis_time, s11a, 'o');
            plot(x_axis_time, s.s11a);
            hold on;
            plot(x_axis_time, s_filtered.s11a);
            for i = 1:length(dS.s11a)
                text((T(i)), peaks.s11a(i), {num2str(dS.s11a(i)), num2str((T(i)))});
            end
            hold on;
            grid on;
            grid minor;
            % plot(x_axis_time, baseline.s11a);
            plot(x_axis_time, detrended.s11a);
            
            % plot(valleySL.t, valleySL.s11a, 'ko');
            % plot(valleySR.t, valleySR.s11a, 'ko');
            
            
            % plot(x_axis_time, bn.s11a, ' r');
            ylabel('\angleS11 in degrees');
            xlabel('time in seconds');
            
            %% plot s21m with annotations
            subplot(2,2,3);
            % plot(x_axis_time, s21m,'o');
            plot(x_axis_time, s.s21m);
            hold on;
            plot(x_axis_time, s_filtered.s21m);
            for i = 1:length(dS.s21m)
                text((T(i)), peaks.s21m(i), {num2str(dS.s21m(i)), num2str((T(i)))});
            end
            hold on;
            grid on;
            grid minor;
            % plot(x_axis_time, baseline.s21m);
            plot(x_axis_time, detrended.s21m);
            
            % plot(valleySL.t, valleySL.s21m, 'ko');
            % plot(valleySR.t, valleySR.s21m, 'ko');
            
            % plot(x_axis_time, bn.s21m, ' r');
            
            ylabel('|S21| in dB');
            xlabel('time in seconds');
            
            %% plot s21a with annotaions
            subplot(2,2,4);
            % plot(x_axis_time, s21a,'o');
            plot(x_axis_time, s.s21a);
            hold on;
            plot(x_axis_time, s_filtered.s21a);
            for i = 1:length(dS.s21a)
                text((T(i)), peaks.s21a(i), {num2str(dS.s21a(i)), num2str((T(i)))});
            end
            hold on;
            grid on;
            grid minor;
            % plot(x_axis_time, baseline.s21a);
            plot(x_axis_time, detrended.s21a);
            
            % plot(valleySL.t, valleySL.s21a, 'ko');
            % plot(valleySR.t, valleySR.s21a, 'ko');
            
            
            % plot(x_axis_time, bn.s21a, ' r');
            ylabel('\angleS21 in degrees');
            xlabel('time in seconds');
            
            %% lock x axis
            ax1 = subplot(2,2,1);
            ax2 = subplot(2,2,2);
            ax3 = subplot(2,2,3);
            ax4 = subplot(2,2,4);
            
            linkaxes([ax1,ax2,ax3,ax4],'x');
            temp = 'filter_adjusted_Processed_Data_Plot';
            saveas(gcf, strcat(working_folder, filesep, temp,'.fig'),'fig');
            close all;

        catch
            disp(strcat(working_folder, "-->doesnt have processsed data"));
            % return;
        end

    end
end
