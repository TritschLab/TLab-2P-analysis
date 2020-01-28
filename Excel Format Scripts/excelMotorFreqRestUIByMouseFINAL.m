%% Summary:
% 
% This script creates a cell array that can be pasted into excel,
% containing the Ca2+ transient frequency during rest and movement for motor
% trials and organized by cell, acquisition, date and mouse.
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% 'mousePeakFreq' - cell array, where each row is a cell, and the
% columns have the labels Mouse / Date / Acq Num / Mot Rest Freq dSPNs / Mot
% Rest Freq iSPNs.
% 
% Author: Jeffrey March, 2018

%% Main Code

[trials, pathname] = uigetfile('*.mat','MultiSelect','on');

cd(pathname)

totalFiles = 0;
mousePeakFreq = {};
cellNum = [1,1];
prevMouse = '';
prevDate = '';
numCondits = 1;

roi = {};
statsMats = {[],[]};

for trial = 1:length(trials);
    load(trials{trial});
    totalFiles = totalFiles + 1
    
    if isfield(data.dF1, 'motorFreqRest')
        roi = {data.dF1, data.dF2};
        
        if ~(strcmp(data.imageFile,prevMouse))
            mousePeakFreq{cellNum(1),1} = data.imageFile;
        end
        
        if  ~(strcmp(data.date,prevDate))
            mousePeakFreq{cellNum(1),2} = data.date;
        end
        
        mousePeakFreq{cellNum(1),3} = data.acqNum;
        
        for cellType = 1:2;
            
            for cellROI = 1:length(roi{cellType}.roiList)
                mousePeakFreq{cellNum(cellType) + cellROI - 1, 4 + (cellType-1)*numCondits} = roi{cellType}.motorFreqRest(cellROI,1);
            end
            
            cellNum(cellType) = cellNum(cellType) + length(roi{cellType}.roiList);
            
        end
        
        for cellType = 1:2
            cellNum(cellType) = max(cellNum);
        end
        
        prevMouse = data.imageFile;
        prevDate = data.date;
    end
                
end
