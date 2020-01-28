%% Summary:
% 
% This script opens and updates all selected files to contain the correct
% Ca2+ transient frequency at rest and during movement.
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

if ~iscell(trials)
    tempTrials = trials;
    trials = cell(1);
    trials{1} = tempTrials;
end

cd(pathname)
totalFiles = 0;

MINW = 4;
noiseThresh = 0.4;
stdAway = 5;

roi = {};

for trial = 1:length(trials);
	load(trials{trial});                               
	totalFiles = totalFiles + 1
    
    roi{1} = data.dF1;
    roi{2} = data.dF2;
      
    for cellType = 1:2;
        if isempty(roi{cellType}.dF)
            continue
        end

        roi2 = zeros(size(roi{cellType}.dF));
        % Z-scoring Cells
        for nCell = 1:size(roi{cellType}.dF,1)
            roi2(nCell,:) = smooth(roi{cellType}.dF(nCell,:),5)';
        end
        
        roi{cellType}.peakFreqRestMot = zeros(size(roi{cellType}.dF,1),2);
        numSpikesCell = zeros(size(roi{cellType}.dF,1),2);
        
        for isMoving = [0, 1];   
            roiData = data;
            
            if isMoving == 0;
                roiData.indOnsets = data.indOnsetsRest;
                roiData.indOffsets = data.indOffsetsRest;
            end

            totalTime = sum(roiData.indOffsets-roiData.indOnsets)/roiData.framerate;
            
            for onset = 1:length(roiData.indOnsets)                    
                onsetMat = roi2(:,roiData.indOnsets(onset):roiData.indOffsets(onset));                
                for nCell = 1:size(onsetMat,1)
                    warning('OFF');
                    MPP = mean(roi2(nCell,roi2(nCell,:) < noiseThresh)) + stdAway*std(roi2(nCell,roi2(nCell,:) < noiseThresh));
                    MPH = MPP;
                    [tempPKS, tempLOCS] = findpeaks(onsetMat(nCell,:),'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW);
                    numSpikesCell(nCell, isMoving + 1) =  numSpikesCell(nCell, isMoving + 1) + length(tempLOCS);
                end
            end    
            
            roi{cellType}.peakFreqRestMot(:,isMoving + 1) = numSpikesCell(:,isMoving + 1)/totalTime;
        end
        
    end
    
    data.dF1 = roi{1};
    data.dF2 = roi{2};
    
    save(trials{trial},'data')  
                    
end