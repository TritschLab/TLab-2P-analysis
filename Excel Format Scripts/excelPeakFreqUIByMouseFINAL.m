%% Summary:
% 
% This script creates a cell array that can be pasted into excel,
% containing all the mouse's Ca2+ transient stats, organized by cell, 
% acquisition, date and mouse.
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% 'mousePeakFreq' - cell array, where each row is a cell, and the columns 
% have the labels Mouse / Date / Acq Num / dSPN Freq Rest / dSPN Freq Mov /
% dSPN Amp Rest / dSPN Amp Mov / dSPN Rise Rest / dSPN Rise Mov / dSPN
% Decay Rest / dSPN Decay Mov / then the same for iSPNs
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
numCondits = 8;

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
            mousePeakFreq{cellNum(cellType) + cellROI - 1, 4 + (cellType-1)*numCondits} = roi{cellType}.peakFreqRestMot(cellROI,1);
            mousePeakFreq{cellNum(cellType) + cellROI - 1, 5 + (cellType-1)*numCondits} = roi{cellType}.peakFreqRestMot(cellROI,2);
            mousePeakFreq{cellNum(cellType) + cellROI - 1, 6 + (cellType-1)*numCondits} = roi{cellType}.peakHeights(cellROI,1);
            mousePeakFreq{cellNum(cellType) + cellROI - 1, 7 + (cellType-1)*numCondits} = roi{cellType}.peakHeights(cellROI,2);
            mousePeakFreq{cellNum(cellType) + cellROI - 1, 8 + (cellType-1)*numCondits} = roi{cellType}.riseTimes(cellROI,1);
            mousePeakFreq{cellNum(cellType) + cellROI - 1, 9 + (cellType-1)*numCondits} = roi{cellType}.riseTimes(cellROI,2);
            mousePeakFreq{cellNum(cellType) + cellROI - 1, 10 + (cellType-1)*numCondits} = roi{cellType}.decayTimes(cellROI,1);
            mousePeakFreq{cellNum(cellType) + cellROI - 1, 11 + (cellType-1)*numCondits} = roi{cellType}.decayTimes(cellROI,2);
        end
        
        cellNum(cellType) = cellNum(cellType) + length(roi{cellType}.roiList);
                
    end
    
    for cellType = 1:2
        cellNum(cellType) = max(cellNum);
    end
                
    prevMouse = data.imageFile;
    prevDate = data.date;
                
end
