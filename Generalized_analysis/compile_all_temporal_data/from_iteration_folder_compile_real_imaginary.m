%%
% this folder starts from the sample class type folder
% looks into 
    % iteration1
    % iteration2
    % iiteration3
% then goes into each frequency and each set to compile the data together
%%
clc;
clear;

cf = pwd;
cf_list = get_folders_list(cf);
cellNumber = 1;
%%
k=1;
for i=1:length(cf_list)
    % this section goes into the iteration folder
    if ispc
        iteration_folder = strcat(cf, "\",cf_list(i).name);
    elseif ismac
        iteration_folder = strcat(cf, "/",cf_list(i).name);
    end
    iteration_number = get_iteration(iteration_folder);
    % then each iteration folder will have different frequencies
    frequency_list = get_folders_list(iteration_folder); % this will show the folders for each frequency
    for ii=1:length(frequency_list)
        if ispc
            frequency_folder = strcat(iteration_folder, "\", frequency_list(ii).name);
        elseif ismac
            frequency_folder = strcat(iteration_folder, "/", frequency_list(ii).name);
        end
        data = compile_all_sets(frequency_folder);
        frequency = get_frequency(frequency_folder);
        temporal_data(ii).frequency= frequency;
        temp = iteration_number*ones(length(data), 1);
        if i>1
            temporal_data(ii).data =  vertcat(temporal_data(ii).data,horzcat(temp, data));
        else
            temporal_data(ii).data =  horzcat(temp, data);
        end
        % k=k+1;
    end
    
end
%%
% filter the data --> remove all data that is not cellNumber = 1;
filtered_data = [];

for i = 1:size(temporal_data, 2)
    filtered_data(i).frequency =  temporal_data(i).frequency;
    k=1;
    for j = 1:size(temporal_data(i).data, 1)
        if temporal_data(i).data(j,12) == cellNumber;
            filtered_data(i).data(k,:) =  temporal_data(i).data(j,:);
            k=k+1;
        end
    end
end

temporal_data = [];
temporal_data = filtered_data;

%%
filename_final_data = strcat("real_imag_temporal_data", "_", num2str(cellNumber));
save(strcat(filename_final_data,'.mat'),'temporal_data');



%% all functions


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
        processedData = to_cartesian(processedData, detrended, s_filtered, t);
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
