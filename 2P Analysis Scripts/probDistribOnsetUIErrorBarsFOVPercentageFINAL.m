%% Summary:
% 
% This script calculates and plots the population frequency at onset and
% offset, normalized to the size of a mean FOV. To run it, the user must
% select the pre-defined data structures from a given condition, and the 
% result will be a trace with each FOV having one average onset. The user 
% can also change the bin size and the smoothing factor.
% 
% Properties for paper: 
% 
% smoothing = 3 bins; window = 200 ms
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% 'spikePerFrame' - a cell array containing the onset matrix for each
% cell type, in which each row is a FOV.
% 
% 'preMat' - average frequency from -3 to -1, before onset. Row is FOV,
% column is cell type.
% 
% 'postMat' - average frequency from 1 to 3, after onset. Row is FOV,
% column is cell type.
% 
% 'ratioPrePost' - ratio of pre- to post-onset frequency. Row is FOV,
% column is cell type.
% 
% Summed population frequency, norm to mean FOV, offset and onset traces
% 
% Author: Jeffrey March, 2018

%% Getting Files and Calculating Mean FOV Size
clear all;

[trials, path] = selectFiles();

% Initializing
roi = {};
sizeFOVMat = zeros(length(trials),2);

% Calculating avg FOV size
for trial = 1:length(trials)
    load(trials{trial});
    
    roi{1} = data.dF1;
    roi{2} = data.dF2;
    for cellType = 1:2;
        sizeFOVMat(trial,cellType) = roi{cellType}.totalSPNs;
    end
end
meanFOV = mean(sizeFOVMat,1);

%% Performing the Actual Calculations for [offset, onset]

for isOnset = [false, true]
    
totalFiles = 0;
LOCS = {};
timeHist = {[],[]};
spikeHist = {[],[]};
spikePerFrame = {[],[]};

MINW = 4;
noiseThresh = 0.4;
stdAway = 5;

roi = {};
totalCells = zeros(1,2);
totalBouts = 0;

for trial = 1:length(trials)
    load(trials{trial});
    totalFiles = totalFiles + 1
    
    LOCS{1} = [];
    LOCS{2} = [];
    
    roi{1} = data.dF1;
    roi{2} = data.dF2;
    roiData = data;
    
    totalBouts = length(roiData.indOnsets); %number of bouts
    totalCells(1) = totalCells(1) + length(roiData.indOnsets)*size(roi{1}.dF,1); %number of d1 cells*bouts
    totalCells(2) = totalCells(2) + length(roiData.indOnsets)*size(roi{2}.dF,1); %number of d2 cells*bouts
    
    for cellType = 1:2;
        
        roi2 = zeros(size(roi{cellType}.dF));
        % Z-scoring Cells
        for nCell = 1:size(roi{cellType}.dF,1)
            roi2(nCell,:) = smooth(roi{cellType}.dF(nCell,:),5)'; %zscore(roi{cellType}(nCell,:));
        end
        
        timeBefore = ceil(4*roiData.framerate);
        timeAfter = timeBefore;
        onsetMat = [];
        
        if isOnset
            for onset = 1:length(roiData.indOnsets)
                onsetMat = cat(1,onsetMat,roi2(:,roiData.indOnsets(onset) - timeBefore:roiData.indOnsets(onset) + timeAfter));
            end
        else
            for onset = 1:length(roiData.indOnsets)
                onsetMat = cat(1,onsetMat,roi2(:,roiData.indOffsets(onset) - timeBefore:roiData.indOffsets(onset) + timeAfter));
            end
        end
        
        % Finding peak locations
        for nCell = 1:size(onsetMat,1)
            warning('OFF');
            currentCell = roi2(mod(nCell - 1,size(roi2,1)) + 1,:);
            MPP = mean(currentCell(currentCell < noiseThresh)) + stdAway*std(currentCell(currentCell < noiseThresh));
            MPH = MPP;
            [tempPKS, tempLOCS] = findpeaks(onsetMat(nCell,:),'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW);
            LOCS{cellType} = cat(2,LOCS{cellType},(tempLOCS - timeBefore)/roiData.framerate);
        end
        
        window = -4:.2:4; % the bins that will make up the frequency trace at onset
        [timeHist{cellType}(trial,:),velHistBins] = histcounts(data.timeDF,window); % counting the number of frames/bin
        [spikeHist{cellType}(trial,:),spikeHistBins] = histcounts(LOCS{cellType},window); % counting the spikes/bin
        spikeHist{cellType}(trial,:) = data.framerate*spikeHist{cellType}(trial,:); % converting spikes/bin to (spikes/bin)*(frames/second)
        timeHist{cellType}(trial,:) = roi{cellType}.totalSPNs*totalBouts*timeHist{cellType}(trial,:)/meanFOV(cellType);
        % Overall, the frequency is going to be spikes * all bouts * all cells * (frames / second) * (mean FOV / (all cells * all bouts * frames))
        % This simplifies to (mean FOV * spikes) / second
    end


end

% Creating Frequency Matrix for Each Cell
for cellType = 1:2
    timeHist{cellType}(timeHist{cellType} == 0) = nan; 
    spikePerFrame{cellType} = spikeHist{cellType}./timeHist{cellType};
    
    preMat(:,cellType) = nanmean(spikePerFrame{cellType}(:,(window >= -3.5 & window < -0.5)),2);
    postMat(:,cellType) = nanmean(spikePerFrame{cellType}(:,(window > 1 & window <= 3)),2);
    ratioPrePost(:,cellType) = preMat(:,cellType)./postMat(:,cellType);
end

temp = cat(2,preMat,postMat,ratioPrePost);

%% Plotting

figure;
hold on
smoothing = 3;

for cellType = 1:2

    if cellType == 1
%         shadedErrorBar(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(totalCells(cellType)),smoothing)','r',1); % n = number of cells
        shadedErrorBar(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(size(spikePerFrame{cellType},1)),smoothing)','r',1); % n = number of FOVs
    else
%         shadedErrorBar(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(totalCells(cellType)),smoothing)','g',1); % n = number of cells
        shadedErrorBar(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(size(spikePerFrame{cellType},1)),smoothing)','g',1); % n = number of FOVs
    end
    
    title('Summed Frequency of All SPNs At Onset/Offset (Norm to Mean FOV)');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    ylim([0,4.5]);
   
end

end