%% Summary:
% 
% This script calculates and plots the fraction of SPNs active during a 
% movement bout. To run it, the user must select the pre-defined data 
% structures from a given condition, and the result will be a trace with 
% each FOV having one average trace. The user can also change the bin size 
% and the smoothing factor. The N is the number of FOVs.
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
% 'spikeHist' - a cell array containing the onToOff matrix for each cell 
% type, in which each row is a FOV.
% 
% (# of active SPNs)/(total SPNs) onToOff traces
% 
% Author: Jeffrey March, 2018

%% Actual Code

clear all;
[trials, path] = selectFiles();

totalFiles = 0;
% LOCS = {};
window = 0:.2:6;
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
        
        onsetSpikeHist = zeros(length(roiData.indOnsets),length(window)-1);
        for onset = 1:length(roiData.indOnsets)
            
            range = roiData.indOffsets(onset) - roiData.indOnsets(onset);
            
            if range/roiData.framerate > max(window) + 1
                
                onsetMat = roi2(:,roiData.indOnsets(onset)-timeBefore:roiData.indOffsets(onset)+timeAfter);
                cellSpikeHist = zeros(size(onsetMat,1),length(window)-1);
                
                for nCell = 1:size(onsetMat,1)
                    
                    MPP = mean(roi2(nCell,roi2(nCell,:) < noiseThresh)) + stdAway*std(roi2(nCell,roi2(nCell,:) < noiseThresh));
                    MPH = MPP;
                    [tempPKS, tempLOCS] = findpeaks(onsetMat(nCell,:),'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW);
                    tempLOCS = (tempLOCS - timeBefore)/roiData.framerate;
                    tempLOCS = tempLOCS(tempLOCS > (range/roiData.framerate)/2 - max(window)/2 & tempLOCS < (range/roiData.framerate)/2 + max(window)/2) - ((range/roiData.framerate)/2 - max(window)/2);
                    cellSpikeHist(nCell,:) = histcounts(tempLOCS,window) > 0;
                    
                end
                
                onsetSpikeHist(onset,:) = sum(cellSpikeHist,1)/roi{cellType}.totalSPNs;
                
            end
            
        end     
        
        spikeHist{cellType}(trial,:) = mean(onsetSpikeHist,1); 
        
    end
    
end


%%
figure;
hold on
smoothing = 5;
for cellType = 1:2
%     window = 0:.005:.25;
    if cellType == 1
        shadedErrorBar(window(1:size(spikeHist{cellType},2)), smooth(nanmean(spikeHist{cellType},1),smoothing)', smooth(nanstd(spikeHist{cellType},1)/sqrt(size(spikeHist{cellType},1)),smoothing)','r',1);
    else
        shadedErrorBar(window(1:size(spikeHist{cellType},2)), smooth(nanmean(spikeHist{cellType},1),smoothing)', smooth(nanstd(spikeHist{cellType},1)/sqrt(size(spikeHist{cellType},1)),smoothing)','g',1);
    end
    
%     if cellType == 1
%         plot(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),5)','r')
%     else
%         plot(window(1:size(spikePerFrame{cellType},2)), smooth(nanmean(spikePerFrame{cellType},1),5)','g')
%     end
    title('Fraction of Total SPNs Active At Onset');
    xlabel('Time (s)');
    ylabel('Fraction');
    ylim([0,.006]);
end


