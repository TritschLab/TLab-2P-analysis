%% Summary:
% 
% This script calculates and plots the spike probability distribution at 
% onset and offset, normalized per active cell. To run it, the user must 
% select the pre-defined data structures from a given condition, and the 
% result will be a trace with each cell at onset being averaged into one 
% large trace. The user can also change the bin size and the smoothing 
% factor. Note: the N is the number of total active cells across FOVs.
% 
% Properties for paper: 
% 
% smoothing = 3 bins; window = 100 ms
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% 'spikePerFrame' - a cell array containing the onset matrix for each
% cell type
% 
% spikes/frame, norm to number of active cells, offset and onset traces
% 
% Author: Jeffrey March, 2018


clear all;
[trials, path] = selectFiles();

for isOnset = true; %[false, true]
    
totalFiles = 0;
% LOCS = {};
timeHist = {[],[]};
spikeHist = {[],[]};
spikePerFrame = {[],[]};

MINW = 4;
noiseThresh = 0.4;
stdAway = 5;

roi = {};
totalCells = zeros(1,2);
totalBouts = 0;
counter = [1,1];

for trial = 1:length(trials)
    load(trials{trial});
    totalFiles = totalFiles + 1
    
%     LOCS{1} = [];
%     LOCS{2} = [];
    
    roi{1} = data.dF1.dF;
    roi{2} = data.dF2.dF;
    roiData = data;
    
    totalBouts = totalBouts + length(roiData.indOnsets); %number of bouts
    totalCells(1) = totalCells(1) + length(roiData.indOnsets)*size(roi{1},1); %number of d1 cells*bouts
    totalCells(2) = totalCells(2) + length(roiData.indOnsets)*size(roi{2},1); %number of d2 cells*bouts
    
    for cellType = 1:2;
        
        roi2 = zeros(size(roi{cellType}));
        % Z-scoring Cells
        for nCell = 1:size(roi{cellType},1)
            roi2(nCell,:) = smooth(roi{cellType}(nCell,:),5)'; %zscore(roi{cellType}(nCell,:));
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
            LOCS = (tempLOCS - timeBefore)/roiData.framerate;
            window = -4:.1:4;
            [timeHist{cellType}(counter(cellType),:),velHistBins] = histcounts(data.timeDF,window);
            [spikeHist{cellType}(counter(cellType),:),spikeHistBins] = histcounts(LOCS,window);
            timeHist{cellType}(counter(cellType),:) = timeHist{cellType}(counter(cellType),:);
            counter(cellType) = counter(cellType) + 1;
        end
        
        
    end


end

for cellType = 1:2
    timeHist{cellType}(timeHist{cellType} == 0) = nan; 
    spikePerFrame{cellType} = spikeHist{cellType}./timeHist{cellType};  
end

%%
figure;
hold on
smoothing = 3;
for cellType = 1:2

    if cellType == 1
        shadedErrorBar(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(totalCells(cellType)),smoothing)','r',1);
%         shadedErrorBar(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(size(spikePerFrame{cellType},1)),smoothing)','r',1);
    else
        shadedErrorBar(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(totalCells(cellType)),smoothing)','g',1);
%         shadedErrorBar(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(size(spikePerFrame{cellType},1)),smoothing)','g',1);
    end
    
%     if cellType == 1
%         plot(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),5)','r')
%     else
%         plot(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),5)','g')
%     end
    title('Ca2+ Event Probability At Onset');
    xlabel('Time (s)');
    ylabel('Probability');
    ylim([0,.0035]);
end
% [N,X] = hist(LOCS{2},window);
% hold on
% plot(X(1:length(N)),smooth(N/totalCells(2),5)','Color',[0,0.8,0])
% title('Spike Probability Distribution (Norm to # of Cells)');
% xlabel('Time After Onset (s)')
% ylabel('Spike Probability in Frame')
% legend('D1','D2','Location','NorthWest')

end
