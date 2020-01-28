%% Summary:
% 
% This script opens and updates all selected files to contain the correct
% Ca2+ transient widths.
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
MINW = 4;
noiseThresh = 0.4;
stdAway = 5;

roi = {};
totalFiles = 0;
warning('OFF');

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
        
        peakWidths = zeros(size(roi2,1),2);       
        
        for isMoving = [0, 1];
            roiData = data;
            
            if isMoving == 0;
                roiData.indOnsets = data.indOnsetsRest;
                roiData.indOffsets = data.indOffsetsRest;
            end
            
            boutMat = zeros(1,size(roi2,2));
            for bout = 1:length(roiData.indOnsets)
                boutMat(roiData.indOnsets(bout):roiData.indOffsets(bout)) = 1;
            end
            
            for nCell = 1:size(roi2,1)

                MPP = mean(roi2(nCell,roi2(nCell,:) < noiseThresh)) + stdAway*std(roi2(nCell,roi2(nCell,:) < noiseThresh));
                baseline = mean(roi2(nCell,roi2(nCell,:) < noiseThresh));
                MPH = MPP;
                [PKS, LOCS, W] = findpeaks(roi2(nCell,:), 'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW);
                W = W(W & boutMat(LOCS))/roiData.framerate;
               
                if isempty(W)
                    peakWidths(nCell, isMoving + 1) = nan;
                else
                    peakWidths(nCell, isMoving + 1) = nanmean(W);
                end % if else
                
            end %for nCell
        end %for movement or rest

        roi{cellType}.peakWidths = peakWidths;
                              
    end %for cellType
    
    data.dF1 = roi{1};
    data.dF2 = roi{2};
    
    save(trials{trial},'data')  

end %for trial