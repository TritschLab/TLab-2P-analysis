%% Summary:
% 
% This script deletes the chosen D1 or D2 neurons from the data structure.
% It accomplishes this by removing the fields and then rerunning the parts
% of martaCleanCodeUIFinal that complete the full data analysis. The script
% acts on the open data file and resaves the file in the folder from which
% you run the code.
% 
% Inputs:
% 
% 'toDelete{1}' - array indicating the dSPNs to delete
% 
% 'toDelete{2}' - array indicating the iSPNs to delete
%
% Outputs:
% 
% Updated data structure, re-saved under its original name in the .mat file
% 
% Author: Jeffrey March, 2018

%% Getting ROIs

toDelete = {[],[]};
toDelete{1} = [2];
toDelete{2} = [];

roi = {data.dF1,data.dF2};
for cellType = 1:2
    if ~isempty(toDelete{cellType})
        for i = toDelete{cellType}(fliplr(1:length(toDelete{cellType})))
            roi{cellType}.dF(i,:) = [];
            roi{cellType}.roiList(i) = [];
        end
    end
end

data.dF1 = roi{1};
data.dF2 = roi{2};
clearvars -except data

%% Creating Movement Onset/Offset Data

for condition = 1:3; % Run through all conditions

    switch condition
        case 1
            roi = data.dF1.dF;
        case 2
            roi = data.dF2.dF;
        case 3
            roi = data.dF3.dF;
    end

    %frame information (both behavior and imaging)
    framerate = data.framerate;
    vel = data.vel;
    fRatio = size(roi,2)/length(vel);

    timeThreshold = 4*framerate; %CHANGE number 4,5,6,7.. a seconda dimensione finestra pre post onset
    velThreshold = 0.004; %CHANGE was .004
    minRestTime = 4*framerate; %CHANGE number 4,5,6,7.. a seconda dimensione finestra pre post onset
    minRunTime = 4*framerate; %CHANGE
    behavior = true; %CHANGE
    % signal = vel; %This can be either absolute value or not
    signal = abs(vel); %This can be either absolute value or not

    % Determining onset/offset indices for mirror
    [onsetsBeh, offsetsBeh] = ...
        getOnsetOffset(signal, velThreshold, minRestTime, minRunTime, behavior);   %%Absolute value of vel CHANGE MARTA

    % Making on/offsets adhere to time/length constraints
    offsetsBeh = offsetsBeh(offsetsBeh < length(signal) - timeThreshold); % Making sure last offsets is at least timeThreshold from the end
    onsetsBeh = onsetsBeh(1:length(offsetsBeh)); % Removing onsets that correspond to removed offsets
    offsetsFinal = offsetsBeh((offsetsBeh - onsetsBeh) > timeThreshold); % Making sure onset to offset is at least timeThreshold in length
    onsetsFinal = onsetsBeh((offsetsBeh - onsetsBeh) > timeThreshold); % Making sure onset to offset is at least timeThreshold in length
    offsetsFinal = offsetsFinal(onsetsFinal > timeThreshold); % Making sure first onset is at least timeThreshold from the beginning (corresponding offset)
    onsetsFinal = onsetsFinal(onsetsFinal > timeThreshold); % Making sure first onset is at least timeThreshold from the beginning

    threshold = 2*std(signal(offsetsFinal(1)+round(framerate):offsetsFinal(1)+3*round(framerate))); % Velocity threshold is 2*stdev of noise
    onsetsFinal = iterToMin(signal,onsetsFinal,threshold,true); % Final onsets array
    offsetsFinal = iterToMin(signal,offsetsFinal,threshold,false); % Final offsets array

    % Necessary again, in case there are problems after iterToMin
    offsetsBeh = offsetsBeh(offsetsBeh < length(signal) - timeThreshold); % Making sure last offsets is at least timeThreshold from the end
    onsetsBeh = onsetsBeh(1:length(offsetsBeh)); % Removing onsets that correspond to removed offsets
    offsetsFinal = offsetsBeh((offsetsBeh - onsetsBeh) > timeThreshold); % Making sure onset to offset is at least timeThreshold in length
    onsetsFinal = onsetsBeh((offsetsBeh - onsetsBeh) > timeThreshold); % Making sure onset to offset is at least timeThreshold in length
    offsetsFinal = offsetsFinal(onsetsFinal > timeThreshold); % Making sure first onset is at least timeThreshold from the beginning (corresponding offset)
    onsetsFinal = onsetsFinal(onsetsFinal > timeThreshold); % Making sure first onset is at least timeThreshold from the beginning

    timeBefore = 4*framerate; % Time before onset/offset (coeff is in seconds)
    timeAfter = 4*framerate; % Time after onset/offset (coeff is in seconds)
    dffOnsets = round(onsetsFinal*fRatio);
    dffOffsets = round(offsetsFinal*fRatio);

    avgdF = mean(roi,1); %Avg dF for roi (roi = dF1 or roi = dF2)
    onsetsMatrixBeh = []; %Initializing
    offsetsMatrixBeh = []; %Initializing
    onsetsMatrixDF = []; %Initializing
    offsetsMatrixDF = []; %Initializing
    onsetToOffsetDF = {}; %Initializing
    onsetToOffsetBeh = {}; %Initializing
    onsetToOffsetTime = {}; %Initializing

    for n = 1:length(dffOnsets)
        onsetsMatrixDF(n,:) = avgdF(dffOnsets(n) - ceil(timeBefore*fRatio):dffOnsets(n) + ceil(timeAfter*fRatio));
        offsetsMatrixDF(n,:) = avgdF(dffOffsets(n) - ceil(timeBefore*fRatio):dffOffsets(n) + ceil(timeAfter*fRatio));
        onsetsMatrixBeh(n,:) = signal(onsetsFinal(n) - ceil(timeBefore):onsetsFinal(n) + ceil(timeAfter));
        offsetsMatrixBeh(n,:) = signal(offsetsFinal(n) - ceil(timeBefore):offsetsFinal(n) + ceil(timeAfter));
        onsetToOffsetDF{n} = avgdF(dffOnsets(n) - ceil(timeBefore*fRatio):dffOffsets(n) + ceil(timeAfter*fRatio));
        onsetToOffsetBeh{n} = signal(onsetsFinal(n) - ceil(timeBefore):offsetsFinal(n) + ceil(timeAfter));
        onsetToOffsetTime{n} = (-ceil(timeBefore*fRatio):length(onsetToOffsetDF{n}) - 1 - ceil(timeAfter*fRatio))/framerate;
    end

    data.indOnsets = onsetsFinal; % Final onset indices
    data.indOffsets = offsetsFinal; % Final offset indices
    data.numBouts = length(onsetsFinal); % Number of bouts
    data.avgBoutDuration = mean(offsetsFinal - onsetsFinal)/framerate; % Avg bout duration
    data.stdBoutDuration = std(offsetsFinal - onsetsFinal)/framerate; % STD bout duration
    data.timeDF = (-ceil(timeBefore*fRatio):ceil(timeAfter*fRatio))/framerate;

    % Saving data into structure
    switch condition
        case 1
            tempstruct = data.dF1;
        case 2
            tempstruct = data.dF2;
        case 3
            tempstruct = data.dF3;
    end


    tempstruct.onsetsMatrixDF = onsetsMatrixDF;
    tempstruct.offsetsMatrixDF = offsetsMatrixDF;
    tempstruct.onsetToOffsetDF = onsetToOffsetDF;
    tempstruct.onsetsMatrixBeh = onsetsMatrixBeh;
    tempstruct.offsetsMatrixBeh = offsetsMatrixBeh;
    tempstruct.onsetToOffsetBeh = onsetToOffsetBeh;
    tempstruct.onsetToOffsetTime = onsetToOffsetTime;


    switch condition
        case 1
            data.dF1 = tempstruct;
        case 2
            data.dF2 = tempstruct;
        case 3
            data.dF3 = tempstruct;
    end

end

clearvars -except data

%% Getting rest onset and offset

framerate = data.framerate;
timeThreshold = 4*framerate; %CHANGE number 4,5,6,7.. a seconda dimensione finestra pre post onset
velThreshold = 0.002; %CHANGE was .004
minRestTime = 4*framerate; %CHANGE number 4,5,6,7.. a seconda dimensione finestra pre post onset
minRunTime =   1; %CHANGE
behavior = true; %CHANGE
signal = abs(data.vel); %This can be either absolute value or not
timeShift = round(.5*framerate);

[onsetsBeh, offsetsBeh] = getOnsetOffset(-signal, -velThreshold, minRunTime, minRestTime, behavior);   %%MARTA CHANGES

% Making on/offsets adhere to time/length constraints
offsetsBeh = offsetsBeh(offsetsBeh < length(signal) - timeThreshold); % Making sure last offsets is at least timeThreshold from the end
onsetsBeh = onsetsBeh(1:length(offsetsBeh)); % Removing onsets that correspond to removed offsets
offsetsFinal = offsetsBeh((offsetsBeh - onsetsBeh) > timeThreshold); % Making sure onset to offset is at least timeThreshold in length
onsetsFinal = onsetsBeh((offsetsBeh - onsetsBeh) > timeThreshold); % Making sure onset to offset is at least timeThreshold in length
offsetsFinal = offsetsFinal(onsetsFinal > timeThreshold); % Making sure first onset is at least timeThreshold from the beginning (corresponding offset)
onsetsFinal = onsetsFinal(onsetsFinal > timeThreshold); % Making sure first onset is at least timeThreshold from the beginning

onsetsFinal = onsetsFinal + timeShift;
offsetsFinal = offsetsFinal - timeShift;

data.indOnsetsRest = onsetsFinal;
data.indOffsetsRest = offsetsFinal;


clearvars -except data

%% Calculating peak frequency at rest and at motion

MINW = 4;
noiseThresh = 0.4;
stdAway = 5;

roi = {};
roi{1} = data.dF1;
roi{2} = data.dF2;

for cellType = 1:2;
    if isempty(roi{cellType}.dF)
        continue
    end
    
    roi2 = zeros(size(roi{cellType}.dF));
    % Z-scoring Cells
    for nCell = 1:size(roi{cellType}.dF,1)
        roi2(nCell,:) = smooth(roi{cellType}.dF(nCell,:),5)';
    end
    
    roi{cellType}.peakFreqRestMot = zeros(size(roi{cellType}.dF,1),2);
    numSpikesCell = zeros(size(roi{cellType}.dF,1),2);
    
    for isMoving = [0, 1];
        roiData = data;
        
        if isMoving == 0;
            roiData.indOnsets = data.indOnsetsRest;
            roiData.indOffsets = data.indOffsetsRest;
        end
        
        totalTime = sum(roiData.indOffsets-roiData.indOnsets)/roiData.framerate;
        
        for onset = 1:length(roiData.indOnsets)
            onsetMat = roi2(:,roiData.indOnsets(onset):roiData.indOffsets(onset));
            for nCell = 1:size(onsetMat,1)
                warning('OFF');
                MPP = mean(roi2(nCell,roi2(nCell,:) < noiseThresh)) + stdAway*std(roi2(nCell,roi2(nCell,:) < noiseThresh));
                MPH = MPP;
                [tempPKS, tempLOCS] = findpeaks(onsetMat(nCell,:),'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW);
                numSpikesCell(nCell, isMoving + 1) =  numSpikesCell(nCell, isMoving + 1) + length(tempLOCS);
            end
        end
        
        roi{cellType}.peakFreqRestMot(:,isMoving + 1) = numSpikesCell(:,isMoving + 1)/totalTime;
    end
    
end

data.dF1 = roi{1};
data.dF2 = roi{2};
    
clearvars -except data

%% Calculating peak height, rise time, decay time, and peak width

MINW = 4;
noiseThresh = 0.4;
stdAway = 5;

PKS = {};
LOCS = {};
roi = {};
warning('OFF');
                
roi{1} = data.dF1;
roi{2} = data.dF2;

for cellType = 1:2;
    if isempty(roi{cellType}.dF)
        continue
    end
    
    roi2 = zeros(size(roi{cellType}.dF));
    % Z-scoring Cells
    for nCell = 1:size(roi{cellType}.dF,1)
        roi2(nCell,:) = smooth(roi{cellType}.dF(nCell,:),5)';
    end
    
    peakHeights = zeros(size(roi2,1),2);
    riseTimes = zeros(size(roi2,1),2);
    decayTimes = zeros(size(roi2,1),2);
    
    for isMoving = [0, 1];
        roiData = data;
        if isMoving == 0;
            roiData.indOnsets = data.indOnsetsRest;
            roiData.indOffsets = data.indOffsetsRest;
        end
        
        boutMat = zeros(1,size(roi2,2));
        for bout = 1:length(roiData.indOnsets)
            boutMat(roiData.indOnsets(bout):roiData.indOffsets(bout)) = 1;
        end
        
        for nCell = 1:size(roi2,1)
            
            MPP = mean(roi2(nCell,roi2(nCell,:) < noiseThresh)) + stdAway*std(roi2(nCell,roi2(nCell,:) < noiseThresh));
            baseline = mean(roi2(nCell,roi2(nCell,:) < noiseThresh));
            MPH = MPP;
            [PKS, LOCS, W] = findpeaks(roi2(nCell,:), 'MinPeakProminence', MPP, 'MinPeakHeight', MPH, 'MinPeakWidth', MINW);
            
            W = W(W & boutMat(LOCS))/roiData.framerate;
            PKS = PKS(PKS & boutMat(LOCS));
            LOCS = LOCS(LOCS & boutMat(LOCS));
            
            decayThreshArray = PKS/3;
            peak80 = PKS*.8;
            peak20 = PKS*.2;
            PKS = roi{cellType}.dF(nCell,LOCS);
            
            if isempty(PKS)
                peakWidths(nCell, isMoving + 1) = nan;
                peakHeights(nCell, isMoving + 1) = nan;
                riseTimes(nCell, isMoving + 1) = nan;
                decayTimes(nCell, isMoving + 1) = nan;
            else
                peakWidths(nCell, isMoving + 1) = nanmean(W);
                peakHeights(nCell,isMoving + 1) = nanmean(PKS);
                riseTimes(nCell,isMoving + 1) = nanmean((iterToMin(roi2(nCell,:),LOCS,peak80,true) - iterToMin(roi2(nCell,:),LOCS,peak20,true))/roiData.framerate);
                decayTimes(nCell,isMoving + 1) = nanmean((iterToMin(roi2(nCell,:),LOCS,decayThreshArray,false) - LOCS)/roiData.framerate);
            end % if else
            
        end %for nCell
    end %for movement or rest
    
    roi{cellType}.peakHeights = peakHeights;
    roi{cellType}.riseTimes = riseTimes;
    roi{cellType}.decayTimes = decayTimes;
    roi{cellType}.peakWidths = peakWidths;
    
    
end %for cellType

data.dF1 = roi{1};
data.dF2 = roi{2};

clearvars -except data

%% Average DF While Moving/Rest or in Bins

for condition = 1:3;
    
    switch condition
        case 1
            roi = data.dF1;
        case 2
            roi = data.dF2;
        case 3
            roi = data.dF3;
    end

    speed = abs(data.vel);
    threshold = .001;
    bins = [0,0.001:0.005:1.001];

    % Calculate binarized data: moving or resting
    moving = speed >= threshold;
    resting = speed < threshold;

    % Calculate mean DF/F when moving or resting
    roi.meanMovDF = nanmean((roi.dF*moving')/sum(moving));
    roi.meanRestDF = nanmean((roi.dF*resting')/sum(resting));
    roi.meanTotDF = nanmean(mean(roi.dF));

    % Calculate mean speed when moving or resting
    roi.meanMovVel = (speed*moving')/sum(moving);
    roi.meanRestVel = (speed*resting')/sum(resting);
    roi.meanTotVel = mean(speed);
    
    % Calculate binarized data: moving or resting
    roi.meanBinDF =[];
    roi.meanBinVel =[];   
    for i = 1:length(bins)-1;
        binNum{i} = speed < bins(i+1) & speed >= bins(i);
        % Calculate mean DF/F when moving or resting
        roi.meanBinDF(i) = nanmean((roi.dF*binNum{i}')/sum(binNum{i}));
        % Calculate mean speed when moving or resting
        roi.meanBinVel(i) = (speed*binNum{i}')/sum(binNum{i});
    end
    
    switch condition
        case 1
            data.dF1 = roi;
        case 2
            data.dF2 = roi;
        case 3
            data.dF3 = roi;
    end
    
end        

clearvars -except data
%% Saving Files

% cd('C:\MATLAB\Calcium Data\Lesion\sChronic regions\')
save([data.imageFile,'_',data.date,'_Acq',num2str(data.acqNum),'.mat'],'data')

clear all