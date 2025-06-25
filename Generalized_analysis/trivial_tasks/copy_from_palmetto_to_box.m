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


% this starts from the folder where you have all the frequencies

clc;
clear;
% sf --> source folder
% df --> destination folder
sf = '/mnt/project/rfcells/December/16th December AR Drug Test CHO';
df = '/media/mdsaifi/Saiful_Lab/CHO_Drug_Effect/December/16th December AR Drug Test CHO';

class_type = 'iteration1';

sf = fullfile(sf,class_type);
df = fullfile(df,class_type);

file1 = "processed_data.mat";
% file2 = "filter_adjusted_Processed_Data_Plot.fig";
file2 = "Processed_Data_Plot.fig";
file3 = "detrend_data.mat";


sf_list = get_folders_list(sf);

for i=1:length(sf_list)

        working_folder = strcat(sf, filesep, sf_list(i).name);
        cd_working_folder = strcat(df, filesep, sf_list(i).name);


        wf_list = get_folders_list(working_folder);
        for ii=1:length(wf_list)
            
            wf = strcat(working_folder, filesep, wf_list(ii).name, filesep);
            cf = strcat(cd_working_folder, filesep, wf_list(ii).name,filesep);
    
            status = copyfile(strcat(wf, file1), cf);
            status = copyfile(strcat(wf, file2), cf);
            % status = copyfile(strcat(wf, file3), cf);
            clearvars wf cf
        end
    % end
end