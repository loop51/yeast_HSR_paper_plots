%%
%  this one compiles from the iteration folders
clc;
clear;

cf = pwd;
cf_list = get_folders_list(cf);
cellNumber = 1;
iteration_number =1;
%%
k=1;

frequency_list = get_folders_list(cf); % this will show the folders for each frequency
for ii=1:length(frequency_list)

    frequency_folder = strcat(cf, filesep, frequency_list(ii).name);

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
    
%%
% filter the data --> remove all data that is not cellNumber = 1;
filtered_data = [];

for i = 1:size(temporal_data, 2)
    filtered_data(i).frequency =  temporal_data(i).frequency;
    k=1;
    for j = 1:size(temporal_data(i).data, 1)
        if temporal_data(i).data(j,12) > 0;
            filtered_data(i).data(k,:) =  temporal_data(i).data(j,:);
            k=k+1;
        end
    end
end

temporal_data = [];
temporal_data = filtered_data;

%%
filename_final_data = strcat("mag_phase_temporal_data", "_", num2str(cellNumber));
save(strcat(filename_final_data,'.mat'),'temporal_data');


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


%%
% all functions

function [data] = compile_all_sets(dl)
    dl_list = get_folders_list(dl);
    data = [];
    for i=1:length(dl_list)
        if ispc
            working_folder = strcat(dl, "\", dl_list(i).name);
            try
            load(strcat(working_folder, "\", "processed_data.mat"));
            data = [data; processedData];
            catch
                disp(strcat(working_folder, "-->doesnt have processsed data"));
            end
        elseif ismac || isunix
            working_folder = strcat(dl, "/", dl_list(i).name);
            try
            load(strcat(working_folder, "/", "processed_data.mat"));
            data = [data; processedData];
            catch
                disp(strcat(working_folder, "-->doesnt have processsed data"));
            end
        end

        

    end
end
%%
function [f] = get_frequency(frequency_directory)
    x = frequency_directory;
    if ismac || isunix
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
    if ismac || isunix
        x = regexp(x, '/', 'split');
        filename_final_data = string(strcat(x(length(x))));
    elseif ispc
        x = regexp(x, '\', 'split');
        filename_final_data = string(strcat(x(length(x))));
    end
    i = str2num(regexp(filename_final_data, '\d+(\.\d+)?', 'match'));
end