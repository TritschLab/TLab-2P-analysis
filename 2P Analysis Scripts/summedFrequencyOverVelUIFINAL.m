%% Summary:
% 
% This script calculates and plots the population frequency, normalized to
% the size of a mean FOV, over velocity. To run it, the user must select
% the pre-defined data structures from a given condition, and the result 
% will be a mean trace with each FOV having one velocity trace. The user
% can also change the velocity bin size and the smoothing factor.
% 
% Properties for paper: 
% 
% smoothing = 5 bins; velWindow = 0.04 m/s; binSize = 0.2 s
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% 'fracOverVel' - array in which the entry number corresponds to a velocity
% bin in 'velWindow', and each value is the population frequency,
% normalized to the mean FOV size.
% 
% Summed population frequency over velocity trace, norm to mean FOV size
% 
% Author: Jeffrey March, 2018

%% Calling each file and extracting FOV size

clear all;

[trials, path] = selectFiles();

roi = {};
sizeFOVMat = zeros(length(trials),2);

%Calculating avg FOV size
for trial = 1:length(trials)
    load(trials{trial});
    
    roi{1} = data.dF1;
    roi{2} = data.dF2;
    for cellType = 1:2;
        sizeFOVMat(trial,cellType) = roi{cellType}.totalSPNs;
    end
end

meanFOV = mean(sizeFOVMat,1);

%% Actual summed frequency part

totalFiles = 0;
binSize = 0.2; % time bin

MINW = 4;
noiseThresh = 0.4;
stdAway = 5;
velWindow = 0:.04:.24; % velocity bin
fracOverVel = {zeros(length(trials),length(velWindow)-1),zeros(length(trials),length(velWindow)-1)};
roi = {};

for trial = 1:length(trials)
    load(trials{trial});
    totalFiles = totalFiles + 1
    
    window = round(0:binSize*data.framerate:length(data.behFrame));
    binnedVel = zeros(1,length(window) - 1);
    for bin = 1:length(window) - 1
        binnedVel(:,bin) = mean(abs(data.vel(:,window(bin) + 1:window(bin + 1))),2);
    end
       
    LOCS{1} = [];
    LOCS{2} = [];
    
    roi{1} = data.dF1;
    roi{2} = data.dF2;
    roiData = data;
    
    for cellType = 1:2;
        
        roi2 = zeros(size(roi{cellType}.dF));   
        % Z-scoring Cells
        for nCell = 1:size(roi{cellType}.dF,1)
            roi2(nCell,:) = smooth(roi{cellType}.dF(nCell,:),5)';
        end
        
        % Finding spike histogram for each cell
        cellSpikeHist = zeros(size(roi2,1),length(window)-1);
        for nCell = 1:size(roi2,1)
            MPP = mean(roi2(nCell,roi2(nCell,:) < noiseThresh)) + stdAway*std(roi2(nCell,roi2(nCell,:) < noiseThresh));
            MPH = MPP;
            [tempPKS, tempLOCS] = findpeaks(roi2(nCell,:),'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW);
            cellSpikeHist(nCell,:) = histcounts(tempLOCS,window); % Number of spikes in each bin (columns) for each cell (rows)
        end        
        spikeHist = data.framerate*sum(cellSpikeHist,1); % for each bin: spikes * frames / second
        
        [timeHist,velHistBins] = histcounts(1:length(data.frameTime),window); % Total number of frames in each bin
        timeHist = roi{cellType}.totalSPNs*timeHist/meanFOV(cellType); % frames * totalSPNs / meanFOV
        
        freqHist = spikeHist./timeHist; % Population frequency
        
        % Each frequency in the corresponding velocity bin
        for velBin = 1:length(velWindow)-1
            fracOverVel{cellType}(trial,velBin) = mean(freqHist(binnedVel >= velWindow(velBin) & binnedVel < velWindow(velBin + 1)));
        end
       
    end  
    
end

%% Plotting

figure;
hold on
smoothing = 5;
for cellType = 1:2
    
    if cellType == 1
        shadedErrorBar(velWindow(1:size(fracOverVel{cellType},2)), smooth(nanmean(fracOverVel{cellType},1),smoothing)', smooth(nanstd(fracOverVel{cellType},1)/sqrt(size(fracOverVel{cellType},1)),smoothing)','r',1);
    else
        shadedErrorBar(velWindow(1:size(fracOverVel{cellType},2)), smooth(nanmean(fracOverVel{cellType},1),smoothing)', smooth(nanstd(fracOverVel{cellType},1)/sqrt(size(fracOverVel{cellType},1)),smoothing)','g',1);
    end
    
    title('Summed Frequency of All SPNs (Norm to Mean FOV) Over Velocity');
    xlabel('Velocity (m/s)');
    ylabel('Frequency (Hz)');
    ylim([0,6]);
    
end