%% Summary:
% 
% This script calculates and plots the fraction of SPNs active at onset and 
% offset. To run it, the user must select the pre-defined data structures 
% from a given condition, and the result will be a trace with each FOV 
% having one average onset. The user can also change the bin size and the 
% smoothing factor. The N is the number of FOVs. Window can be changed to be 
% defined by "fraction of cells active during the 2 secs before onset and 
% the 2 secs after onset." Can also be changed between "fraction of total
% SPNs" and "fraction of active SPNs.
% 
% Properties for paper: 
% 
% smoothing = 5 bins; window = 200 ms
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% 'spikeHist' - a cell array containing the onset matrix for each
% cell type, in which each row is a FOV.
% 
% (# of active SPNs)/(total SPNs) offset and onset traces
% 
% Author: Jeffrey March, 2018

%% Actual Code

clear all;
[trials, path] = selectFiles();

for isOnset = [true, false]
    
totalFiles = 0;
warning('off')
% LOCS = {};
window = -4:.2:4; %smooth graph
% window = [-3,-1,1,3]; %pre-post graph
spikeHist = {zeros(length(trials),length(window)-1),zeros(length(trials),length(window)-1)};

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
    
    totalBouts = totalBouts + length(roiData.indOnsets); %number of bouts
    totalCells(1) = totalCells(1) + length(roiData.indOnsets)*size(roi{1},1); %number of d1 cells*bouts
    totalCells(2) = totalCells(2) + length(roiData.indOnsets)*size(roi{2},1); %number of d2 cells*bouts
    
    for cellType = 1:2;
        
        roi2 = zeros(size(roi{cellType}.dF));
        % Z-scoring Cells
        for nCell = 1:size(roi{cellType}.dF,1)
            roi2(nCell,:) = smooth(roi{cellType}.dF(nCell,:),5)';
        end
        
        timeBefore = ceil(4*roiData.framerate);
        timeAfter = timeBefore;
        onsetMat = [];
        
        if isOnset
            onsets = roiData.indOnsets;
        else
            onsets = roiData.indOffsets;
        end
        
        onsetSpikeHist = zeros(length(onsets),length(window)-1);
        for onset = 1:length(onsets)          
            onsetMat = roi2(:,onsets(onset)-timeBefore:onsets(onset)+timeAfter);
            cellSpikeHist = zeros(size(onsetMat,1),length(window)-1);
            for nCell = 1:size(onsetMat,1)
                MPP = mean(roi2(nCell,roi2(nCell,:) < noiseThresh)) + stdAway*std(roi2(nCell,roi2(nCell,:) < noiseThresh));
                MPH = MPP;
                [tempPKS, tempLOCS] = findpeaks(onsetMat(nCell,:),'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW);
                tempLOCS = (tempLOCS - timeBefore)/roiData.framerate;
                cellSpikeHist(nCell,:) = histcounts(tempLOCS,window) > 0;
            end
%             onsetSpikeHist(onset,:) = sum(cellSpikeHist,1)/roi{cellType}.totalSPNs; % For fraction of total neurons
            onsetSpikeHist(onset,:) = sum(cellSpikeHist,1)/size(roi{cellType}.dF,1); % For fraction of active neurons
        end      
        spikeHist{cellType}(trial,:) = mean(onsetSpikeHist,1);
        
%         preMat(:,cellType) = nanmean(spikeHist{cellType}(:,(window >= -3.5 & window < -0.5)),2);
%         postMat(:,cellType) = nanmean(spikeHist{cellType}(:,(window > 1 & window <= 3)),2);
%         ratioPrePost(:,cellType) = preMat(:,cellType)./postMat(:,cellType);
%         
    end
end

% temp = cat(2,preMat,postMat,ratioPrePost);

%% Plotting

figure;
hold on
smoothing = 5; %smoothing = 5 for continuous graph; 1 for not continuous
for cellType = 1:2

    if cellType == 1
        shadedErrorBar(window(1:size(spikeHist{cellType},2)), smooth(nanmean(spikeHist{cellType},1),smoothing)', smooth(nanstd(spikeHist{cellType},1)/sqrt(size(spikeHist{cellType},1)),smoothing)','r',1);
    else
        shadedErrorBar(window(1:size(spikeHist{cellType},2)), smooth(nanmean(spikeHist{cellType},1),smoothing)', smooth(nanstd(spikeHist{cellType},1)/sqrt(size(spikeHist{cellType},1)),smoothing)','g',1);
    end
    
    title('Fraction of Total SPNs Active At Onset');
    xlabel('Time (s)');
    ylabel('Fraction');
%     ylim([0,.05]); % For fraction of total SPNs, pre/post
%     ylim([0,.006]); % For fraction of total SPNs, smooth curve
    
    
end

end
