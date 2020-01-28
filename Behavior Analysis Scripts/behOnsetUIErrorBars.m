%% Summary:
% 
% This script calculates and plots average speed traces at onset, from
% onset to offset, and at offset. The N is the number of FOVs.
% 
% Properties for paper: 
% 
% window = 6 s
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% 'onsetMat' - a cell array containing the onset matrix for each
% cell type, in which each row is a FOV.
% 
% 'onToOffMat' - average frequency from -3 to -1, before onset. Row is FOV,
% column is cell type.
% 
% 'offsetMat' - average frequency from 1 to 3, after onset. Row is FOV,
% column is cell type.
% 
% average speed traces at onset, onToOff, and offset
% 
% Author: Jeffrey March, 2018

%% Creating onset, OnToOff, and Offset Matrices

[trials, path] = selectFiles();

totalFiles = 0;
onsetMat = [];
onToOffMat = [];
offsetMat = [];

roi = {};
window = 6;

for trial = 1:length(trials)
    load(trials{trial});
    totalFiles = totalFiles + 1
    
    roi{1} = data.dF1;
    roi{2} = data.dF2;
    roiData = data;
    
    timeBefore = ceil((window/2)*roiData.framerate);
    timeAfter = timeBefore;
    
    cellType = 1;
        
    onToOff = [];
    for onset = 1:length(roiData.indOnsets)
        range = roiData.indOffsets(onset) - roiData.indOnsets(onset);
        if range/roiData.framerate > window + 1
            midpoint = round((roiData.indOnsets(onset)+roiData.indOffsets(onset))/2);
            onToOff = cat(1,onToOff,mean(abs(roiData.vel(:,(midpoint-timeBefore):(midpoint + timeBefore))),1));
        end
    end
    
    onsetMat = cat(1,onsetMat,mean(roi{1}.onsetsMatrixBeh,1));
    onToOffMat = cat(1,onToOffMat,mean(onToOff,1));
    offsetMat = cat(1,offsetMat,mean(roi{1}.offsetsMatrixBeh,1));
    
    
end

%% Plotting Onsets

figure;
shadedErrorBar(data.timeDF, nanmean(onsetMat,1), nanstd(onsetMat,1)/sqrt(size(onsetMat,1)),'k',1);
title('Average Speed At Onset');
xlabel('Time (s)');
ylabel('Speed (m/s)');
ylim([0,.10])

%% Plotting On To Off

figure;
timeOnToOff = 1/ceil(data.framerate):1/ceil(data.framerate):size(onToOff,2)/ceil(data.framerate);
shadedErrorBar(timeOnToOff, nanmean(onToOffMat,1), nanstd(onToOffMat,1)/sqrt(size(onToOffMat,1)),'k',1);
title('Average Speed During Bout');
xlabel('Time (s)');
ylabel('Speed (m/s)');
ylim([0,.10])


%% Plotting Offsets

figure;
shadedErrorBar(data.timeDF, nanmean(offsetMat,1), nanstd(offsetMat,1)/sqrt(size(offsetMat,1)),'k',1);
title('Average Speed At Offset');
xlabel('Time (s)');
ylabel('Speed (m/s)');
ylim([0,.10])


