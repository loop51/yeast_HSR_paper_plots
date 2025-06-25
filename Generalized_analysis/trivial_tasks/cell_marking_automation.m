clear;
close all;

load('processed_data.mat');

openfig("filter_adjusted_Processed_Data_Plot.fig");

set_param = 1;



for i=1:size(processedData,1)
    processedData(i, 11) = set_param;
    if processedData(i, 5) == 0 
        processedData(i, 11) = 0;
    end
end

save('processed_data.mat','processedData');
clear;

load('processed_data.mat');