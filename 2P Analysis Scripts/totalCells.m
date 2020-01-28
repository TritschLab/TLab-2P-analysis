%% Summary:
% 
% This script counts the total number of cells of each cell type that are 
% analyzed for a given condition.
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% 'totalCells' - a cell array containing the total number of cells of each
% cell type that are analyzed for a given condition
% 
% Author: Jeffrey March, 2018

%% Actual code

[trials, pathname] = selectFiles();

totalCells = {0,0};
totalFiles = 0;
roi = {};

for trial = 1:length(trials)
    load(trials{trial});
    totalFiles = totalFiles + 1
    
    roi{1} = data.dF1.dF;
    roi{2} = data.dF2.dF;
    
    for cellType = 1:2
        totalCells{cellType} = totalCells{cellType} + size(roi{cellType},1);
    end
end
        