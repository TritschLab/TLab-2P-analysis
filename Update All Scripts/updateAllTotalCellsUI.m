%% Summary:
% 
% This script opens and updates all selected files to contain the correct
% total cell count for the FOV.
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

totalFiles = 0;

% Manually create totalCellMat from 'Total D1s and D2s Across Conditions'
% Excel sheet

for trial = 1:length(trials);
    if length(trials) ~= size(totalCellMat,1)
        display('Error: not equal lengths');
        break
    end
	load(trials{trial});
    totalFiles = totalFiles + 1
                
    roi{1} = data.dF1;
    roi{2} = data.dF2;
    
    for cellType = 1:2;
        roi{cellType}.totalSPNs = totalCellMat(trial,cellType);
    end
    
    data.dF1 = roi{1};
    data.dF2 = roi{2};
    
    save(trials{trial},'data');
end
    

