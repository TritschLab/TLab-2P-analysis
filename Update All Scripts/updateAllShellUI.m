%% Summary:
% 
% This script contains the shell for all the updateAll* scripts! Just plug
% in what it is that you want to change in the "for cellType" loop.
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

[trials, pathname] = uigetfile('*.mat','MultiSelect','on');

cd(pathname)

if ~iscell(trials)
    tempTrials = trials;
    trials = cell(1);
    trials{1} = tempTrials;
end

totalFiles = 0;

for trial = 1:length(trials);
	load(trials{trial});
    totalFiles = totalFiles + 1
                
    roi{1} = data.dF1;
    roi{2} = data.dF2;
    roiData = data;
                
    for cellType = 1:2;       
                

                              
    end %for cellType
    
    data.dF1 = roi{1};
    data.dF2 = roi{2};
    
    save(trials{trial},'data')  

end %for trial