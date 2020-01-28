%% Summary:
% 
% This script calculates and plots the per cell frequency at onset and
% offset, normalized to the total number of SPNs. To run it, the user must
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
% per cell population frequency, norm to total SPNs, offset and onset traces
% 
% Author: Jeffrey March, 2018

clear all;
[trials, path] = selectFiles();

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
        for nCell = 1:size(onsetMat,1)
            warning('OFF');
            currentCell = roi2(mod(nCell - 1,size(roi2,1)) + 1,:);
            MPP = mean(currentCell(currentCell < noiseThresh)) + stdAway*std(currentCell(currentCell < noiseThresh));
            MPH = MPP;
            [tempPKS, tempLOCS] = findpeaks(onsetMat(nCell,:),'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW);
            LOCS{cellType} = cat(2,LOCS{cellType},(tempLOCS - timeBefore)/roiData.framerate);
        end
        
        window = -4:.2:4;
        [timeHist{cellType}(trial,:),velHistBins] = histcounts(data.timeDF,window);
        [spikeHist{cellType}(trial,:),spikeHistBins] = histcounts(LOCS{cellType},window);
        spikeHist{cellType}(trial,:) = data.framerate*spikeHist{cellType}(trial,:);
        timeHist{cellType}(trial,:) = roi{cellType}.totalSPNs*totalBouts*timeHist{cellType}(trial,:);
        
    end


end

for cellType = 1:2
    timeHist{cellType}(timeHist{cellType} == 0) = nan; 
    spikePerFrame{cellType} = spikeHist{cellType}./timeHist{cellType};
    
    preMat(:,cellType) = nanmean(spikePerFrame{cellType}(:,(window >= -3.5 & window < -0.5)),2);
    postMat(:,cellType) = nanmean(spikePerFrame{cellType}(:,(window > 1 & window <= 3)),2);
    ratioPrePost(:,cellType) = preMat(:,cellType)./postMat(:,cellType);
end

temp = cat(2,preMat,postMat,ratioPrePost);

%%
figure;
hold on
smoothing = 3;
for cellType = 1:2
%     window = 0:.005:.25;
    if cellType == 1
%         shadedErrorBar(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(totalCells(cellType)),smoothing)','r',1);
        shadedErrorBar(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(size(spikePerFrame{cellType},1)),smoothing)','r',1);
    else
%         shadedErrorBar(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(totalCells(cellType)),smoothing)','g',1);
        shadedErrorBar(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(size(spikePerFrame{cellType},1)),smoothing)','g',1);
    end
    
%     if cellType == 1
%         plot(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),5)','r')
%     else
%         plot(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),5)','g')
%     end
    title('Ca2+ Event Probability At Onset');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    ylim([0,.03]);
  
    
end
% [N,X] = hist(LOCS{2},window);
% hold on
% plot(X(1:length(N)),smooth(N/totalCells(2),5)','Color',[0,0.8,0])
% title('Spike Probability Distribution (Norm to # of Cells)');
% xlabel('Time After Onset (s)')
% ylabel('Spike Probability in Frame')
% legend('D1','D2','Location','NorthWest')

end
%%