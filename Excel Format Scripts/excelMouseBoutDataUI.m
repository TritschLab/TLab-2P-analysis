%% Summary:
% 
% This script creates a cell array that can be pasted into excel,
% containing the mouse's movement bout information, organized by bout, 
% acquisition, date and mouse.
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% 'mouseBoutData' - cell array, where each row is a bout, and the columns 
% have the labels Mouse / Date / Acq Num / Bout Duration / Mean Bout 
% Velocity / Peak Bout Velocity
% 
% Author: Jeffrey March, 2018

%% Main Code

[trials, path] = uigetfile('*.mat','MultiSelect','on');
cd(path)

if ~iscell(trials)
    tempTrials = trials;
    trials = cell(1);
    trials{1} = tempTrials;
end

totalFiles = 0;
mouseBoutData = {};
boutNum = 1;
prevMouse = '';
prevDate = '';

for trial = 1:length(trials);
    load(trials{trial});
    totalFiles = totalFiles + 1
                
    if ~(strcmp(data.imageFile,prevMouse))
        mouseBoutData{boutNum,1} = data.imageFile;
    end 
                
    if  ~(strcmp(data.date,prevDate))
        mouseBoutData{boutNum,2} = data.date;
    end
                
    mouseBoutData{boutNum,3} = data.acqNum;
                
    if strcmp(data.imageFile(1),'F')
        roi = data;
    end
    if strcmp(data.imageFile(1),'S')
        roi = data;
        roi.sampleRate = data.framerate;
        roi.FP = data.dF1.dF;
        roi.onsetToOffsetBeh = data.dF1.onsetToOffsetBeh;
        roi.onsetToOffsetTime = data.dF1.onsetToOffsetTime;
    end
    
    boutDurations = [];
    meanVels = [];
    peakVels = [];
    
    fRatio = size(roi.FP,2)/length(roi.vel);
    timeBefore = 4*roi.sampleRate; % Time before onset/offset (coeff is in seconds)
    timeAfter = 4*roi.sampleRate; % Time after onset/offset (coeff is in seconds)
    for bout= 1:length(roi.onsetToOffsetTime)
        boutDurations(bout) = roi.onsetToOffsetTime{bout}(end) - 4; % 4 seconds is time after offset recorded
        meanVels(bout) = mean(abs(roi.onsetToOffsetBeh{bout}(round((timeBefore + 1)*fRatio:end - round((timeBefore + 1)*fRatio)))));
        peakVels(bout) = max(abs(roi.onsetToOffsetBeh{bout}(round((timeBefore + 1)*fRatio:end - round((timeBefore + 1)*fRatio)))));
    end
    
    for bout = 1:length(boutDurations)
        mouseBoutData{boutNum + bout - 1, 4} = boutDurations(bout);
        mouseBoutData{boutNum + bout - 1, 5} = meanVels(bout);
        mouseBoutData{boutNum + bout - 1, 6} = peakVels(bout);
    end
    
    boutNum = boutNum + length(boutDurations);
    
    prevMouse = data.imageFile;
    prevDate = data.date;
end

