%% Summary:
% 
% This script creates a cell array that can be pasted into excel,
% containing each FOV's percent of movement-specific cells, normalized to 
% the number of active cells and organized by acquisition, date and mouse.
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% 'mousePeakFreq' - cell array, where each row is an FOV, and the columns 
% have the labels Mouse / Date / Acq Num / % Rest Specific dSPNs / % Non-
% Specific dSPNs / % Movement Specific dSPNs / then the same for iSPNs.
% 
% Author: Jeffrey March, 2018

%% Main Code

[trials, pathname] = selectFiles;

totalFiles = 0;
mousePeakFreq = {};
cellNum = [1,1];
prevMouse = '';
prevDate = '';
numCondits = 3;

roi = {};
statsMats = {[],[]};

for trial = 1:length(trials);
    load(trials{trial});
    totalFiles = totalFiles + 1
    
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
            mousePeakFreq{cellNum(cellType), 4 + (cellType-1)*numCondits} = sum(roi{cellType}.peakFreqRestMot(:,1) > 0 & roi{cellType}.peakFreqRestMot(:,2) == 0)/sum(roi{cellType}.peakFreqRestMot(:,1) > 0 | roi{cellType}.peakFreqRestMot(:,2) > 0);
            mousePeakFreq{cellNum(cellType), 5 + (cellType-1)*numCondits} = sum(roi{cellType}.peakFreqRestMot(:,1) > 0 & roi{cellType}.peakFreqRestMot(:,2) > 0)/sum(roi{cellType}.peakFreqRestMot(:,1) > 0 | roi{cellType}.peakFreqRestMot(:,2) > 0);
            mousePeakFreq{cellNum(cellType), 6 + (cellType-1)*numCondits} = sum(roi{cellType}.peakFreqRestMot(:,1) == 0 & roi{cellType}.peakFreqRestMot(:,2) > 0)/sum(roi{cellType}.peakFreqRestMot(:,1) > 0 | roi{cellType}.peakFreqRestMot(:,2) > 0);    
        end
        
        cellNum(cellType) = cellNum(cellType) + 1;
                
    end
    
    for cellType = 1:2
        cellNum(cellType) = max(cellNum);
    end
                
    prevMouse = data.imageFile;
    prevDate = data.date;
                
end
