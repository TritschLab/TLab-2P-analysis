% Summary:
% 
% This script calculates and plots the spike probability distribution in the
% middle of the bout, normalized per active cell. To run it, the user must 
% select the pre-defined data structures from a given condition, and the 
% result will be a trace with each FOV having one average onset. The user 
% can also change the bin size and the smoothing factor. The N is the
% number of FOVs
% 
% Properties for paper: 
% 
% smoothing = 3 bins; bins = 200 ms; window = 6s
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
% spikes/frame, norm to number of active cells, mid-bout trace
% 
% Author: Jeffrey March, 2018

%%

clear all;
[trials, path] = uigetfile('*.mat','MultiSelect','on');
cd(path)

totalFiles = 0;
LOCS = {};
timeHist = {[],[]};
spikeHist = {[],[]};
spikePerFrame = {[],[]};

MINW = 4;
noiseThresh = 0.4;
stdAway = 5;

roi = {};
window = 6;
bins = 0:.2:6;

for trial = 1:length(trials)
    load(trials{trial});
    totalFiles = totalFiles + 1
    
    LOCS{1} = [];
    LOCS{2} = [];
    
    roi{1} = data.dF1.dF;
    roi{2} = data.dF2.dF;
    roiData = data;
    
    timeBefore = ceil(4*roiData.framerate);
    timeAfter = timeBefore;
    
    for cellType = 1:2;
        
        roi2 = zeros(size(roi{cellType}));
        % Z-scoring Cells
        for nCell = 1:size(roi{cellType},1)
            roi2(nCell,:) = smooth(roi{cellType}(nCell,:),5)'; %zscore(roi{cellType}(nCell,:));
        end
        
        %                 onsetMat = [];
        for onset = 1:length(roiData.indOnsets)
            range = roiData.indOffsets(onset) - roiData.indOnsets(onset);
            if range/roiData.framerate > window + 1
                onsetMat = roi2(:,roiData.indOnsets(onset)-timeBefore:roiData.indOffsets(onset)+timeAfter);
                for nCell = 1:size(onsetMat,1)
                    MPP = mean(roi2(nCell,roi2(nCell,:) < noiseThresh)) + stdAway*std(roi2(nCell,roi2(nCell,:) < noiseThresh));
                    MPH = MPP;
                    [tempPKS, tempLOCS] = findpeaks(onsetMat(nCell,:),'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW);
                    tempLOCS = (tempLOCS - timeBefore)/roiData.framerate;
                    tempLOCS = tempLOCS(tempLOCS > (range/roiData.framerate)/2 - window/2 & tempLOCS < (range/roiData.framerate)/2 + window/2) - ((range/roiData.framerate)/2 - window/2);
                    LOCS{cellType} = cat(2,LOCS{cellType}, tempLOCS);
                    %                 LOCS{cellType} = cat(2,LOCS{cellType},round(100*tempLOCS/size(onsetMat,2)));
                end            
                
            end
        end
        
        timeHist{cellType}(trial,:) = histcounts((0:(1/data.framerate):6),bins);
        spikeHist{cellType}(trial,:) = histcounts(LOCS{cellType},bins);
        timeHist{cellType}(trial,:) = length(roiData.indOnsets)*size(onsetMat,1)*timeHist{cellType}(trial,:);
        
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
        shadedErrorBar(bins(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(size(spikePerFrame{cellType},1)),smoothing)','r',1);
    else
        shadedErrorBar(bins(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),smoothing)', smooth(nanstd(spikePerFrame{cellType},1)/sqrt(size(spikePerFrame{cellType},1)),smoothing)','g',1);
    end
    
%     if cellType == 1
%         plot(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),5)','r')
%     else
%         plot(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),5)','g')
%     end
    title('Ca2+ Event Probability During Bout');
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


