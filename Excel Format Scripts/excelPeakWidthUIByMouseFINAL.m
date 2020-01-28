%% Summary:
% 
% This script creates a cell array that can be pasted into excel,
% containing each cell's Ca2+ transient width, organized by cell, 
% acquisition, date and mouse.
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% 'mousePeakFreq' - cell array, where each row is a cell, and the columns 
% have the labels Mouse / Date / Acq Num / dSPN Width Rest / dSPN Width Mov
% / iSPN Width Rest / iSPN Width Mov
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
numCondits = 2;

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
            mousePeakFreq{cellNum(cellType) + cellROI - 1, 4 + (cellType-1)*numCondits} = roi{cellType}.peakWidths(cellROI,1);
            mousePeakFreq{cellNum(cellType) + cellROI - 1, 5 + (cellType-1)*numCondits} = roi{cellType}.peakWidths(cellROI,2);

        end
        
        cellNum(cellType) = cellNum(cellType) + length(roi{cellType}.roiList);
                
    end
    
    for cellType = 1:2
        cellNum(cellType) = max(cellNum);
    end
                
    prevMouse = data.imageFile;
    prevDate = data.date;
                
end
