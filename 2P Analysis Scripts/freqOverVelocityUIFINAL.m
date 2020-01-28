%% Summary:
% 
% This script calculates and plots the Ca2+ transient prob distribution per 
% cell, over velocity. To run it, the user must select the pre-defined data
% structures from a given condition, and the result will be a mean trace 
% with each FOV having one velocity trace. The user can also change the 
% velocity bin size and the smoothing factor.
% 
% Properties for paper: 
% 
% smoothing = 5 bins; window = 0.004 m/s
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% 'spikePerFrame' - array in which the entry number corresponds to a velocity
% bin in 'window', and each value is the average spikes per frame per cell.
% 
% average probability distribution per cell over velocity trace.
% 
% Author: Jeffrey March, 2018

%% Code calculations
clear all;

[trials, path] = selectFiles();

totalFiles = 0;
velDistrib = {};
VELS = {};

MINW = 4;
noiseThresh = 0.4;
stdAway = 5;

roi = {};
spikePerFrame = {};

for trial = 1:length(trials);
    load(trials{trial});
    totalFiles = totalFiles + 1
     
    roi{1} = data.dF1.dF;
    roi{2} = data.dF2.dF;
    roiData = data;
    
    for cellType = 1:2;
        
        velDistrib{1} = [];
        velDistrib{2} = [];
        VELS{1} = [];
        VELS{2} = [];
        
        roi2 = zeros(size(roi{cellType}));
        % Z-scoring Cells
        for nCell = 1:size(roi{cellType},1)
            roi2(nCell,:) = smooth(roi{cellType}(nCell,:),5)';
        end
        
        timeBefore = ceil(4*roiData.framerate);
        timeAfter = timeBefore;
        onsetMat = [];
        
        for nCell = 1:size(roi2,1)
            warning('OFF');
            MPP = mean(roi2(nCell,roi2(nCell,:) < noiseThresh)) + stdAway*std(roi2(nCell,roi2(nCell,:) < noiseThresh));
            MPH = MPP;
            [tempPKS, tempLOCS] = findpeaks(roi2(nCell,:),'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW);
            VELS{cellType} = cat(2,VELS{cellType},abs(data.vel(tempLOCS)));
            velDistrib{cellType} = cat(2,velDistrib{cellType},abs(data.vel));
        end
        
        window = 0:.004:.204;
%         window = [0:.001:.004, .006:.002:.01, .02:.01:.25];
        [velHist{cellType}(trial,:),velHistBins] = histcounts(velDistrib{cellType},window);
        [spikeHist{cellType}(trial,:),spikeHistBins] = histcounts(VELS{cellType},window);
    end
    
end

for cellType = 1:2
    velHist{cellType}(velHist{cellType} < 100) = nan; 
    spikePerFrame{cellType} = spikeHist{cellType}./velHist{cellType};  
end

%% Plotting

figure;
hold on

smoothing = 5;
for cellType = 1:2
    
    if cellType == 1
        shadedErrorBar(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(size(spikePerFrame{cellType},1)),smoothing)','r',1);
    else
        shadedErrorBar(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(size(spikePerFrame{cellType},1)),smoothing)','g',1);
    end
    
    title('Ca2+ Event Probability Over Velocity');
    xlabel('Speed (m/s)');
    ylabel('Probability');
    ylim([0,.004]);
    
end