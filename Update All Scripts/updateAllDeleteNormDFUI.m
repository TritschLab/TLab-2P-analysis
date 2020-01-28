%% Summary:
% 
% This script opens and updates all selected files to delete the normDF
% fields.
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% Updated data structure, re-saved under its original name in the .mat file
% 
% Author: Jeffrey March, 2018

%% Main Code

[trials] = uigetfile('*.mat','MultiSelect','on');

cd('C:\MATLAB\Calcium Data\')
totalFiles = 0;

for trial = 1:length(trials);
	load(trials{trial});                               
	totalFiles = totalFiles + 1
    if isfield(data, 'normDF1')
        data = rmfield(data,'normDF1');
        data = rmfield(data,'normDF2');
        data = rmfield(data,'normDF3');
    end
    
    save(trials{trial},'data')  
                    
end