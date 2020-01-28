% Summary:
% 
% This script is used to check that all peaks are getting detected by our 
% findpeaks threshold by plotting every cell's DF/F trace with a stem plot 
% of each peak superimposed on the continuous signal. It operates by
% opening up the data matrix and then running. Note: if there are too many
% cells, you may have run twice, once for dSPNs and once for iSPNs (by
% changing cellType = 1:2 to 1 or 2 individually).
% 
% Properties for paper: 
% 
% MINW = 4 points; noiseThresh = 0.4 DF/F; stdAway = 5 std
% 
% Inputs:
% 
% User-opened data file.
%
% Outputs:
% 
% A figure for every cell's DF/F trace with a stem plot of each peak 
% superimposed on the continuous signal
% 
% Author: Jeffrey March, 2018

%% CORRECT!!

totalFiles = 0;
MINW = 4;
noiseThresh = 0.4; % What is considered noise
stdAway = 5; % How many std from the mean of the noise must a peak be
totalFiles = totalFiles + 1

roi{1} = data.dF1.dF;
roi{2} = data.dF2.dF;
roiData = data;
warning('OFF');

for cellType = 1:2 % Change to run each loop individually if there are too many cells

    roi2 = zeros(size(roi{cellType}));
    % Z-scoring Cells
    for nCell = 1:size(roi{cellType},1)
        roi2(nCell,:) = smooth(roi{cellType}(nCell,:),5)'; %zscore(roi{cellType}(nCell,:));
    end

    for nCell = 1:size(roi2,1)      
        MPP = mean(roi2(nCell,roi2(nCell,:) < noiseThresh)) + stdAway*std(roi2(nCell,roi2(nCell,:) < noiseThresh));
        MPH = MPP;
        [tempPKS, tempLOCS] = findpeaks(roi2(nCell,:),'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW);
        figure;
        plot(roiData.frameTime,roi2(nCell,:))
        hold on;
        stem(roiData.frameTime(tempLOCS),roi2(nCell,tempLOCS))
        title(['Threshold: ', num2str(mean(roi2(nCell,roi2(nCell,:) < noiseThresh)) + stdAway*std(roi2(nCell,roi2(nCell,:) < noiseThresh)))])
    end
    
end
