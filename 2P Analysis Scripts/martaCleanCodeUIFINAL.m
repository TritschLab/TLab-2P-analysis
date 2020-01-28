%% Summary:
% 
% This script calculates and saves all data analyses for 2P analysis. To run,
% you need to define data.imageFile as the mouse to be analyzed, data.date
% to be the date of the acquisition, and firstFile to be the acquisition
% number to be analyzed. You must also make sure the folder with the
% acquisition 2P file, the associated movie files, and the behavior files
% are in the path. Aside from that, this code should run without errors and
% create the final data structure that will be used for all analysis. The
% code is well-commented.
% 
% Properties for paper: 
% 
% wheelRadius = .06 m; velThreshold (mov) = 0.004 m/s;
% velThreshold (rest) = 0.002 m/s; 
% 
% Inputs:
% 
% 'data.imageFile' - mouse name
% 
% 'data.date' - date of acquisition
% 
% 'firstFile' - acquisition number
% 
% acquisition 2P file, related movies, and behavior traces
%
% Outputs:
% 
% 'data' - large structure of data containing all relevant analyses. See
% sheet on the "Code List.xls" spreadsheet entitled "2P Data Struct
% Schematic"
% 
% plots that will not be used
% 
% Author: Jeffrey March, 2018

%% Getting ROIs

clear all;
close all;
clc;

data = struct; % Defining data as an empty structure

data.imageFile = 'SW024'; %CHANGE
data.date = '180502'; %CHANGE

load([data.imageFile,'.mat']);  %Loading the obj file associated with data.imagefile
myObj = eval(data.imageFile); %Bad practice, I know...


[sliceNum, channelNum] = getDF(myObj); %Opening object and prepping for extraction

roiGroup = 1; %Getting D1 Neuron Data

% Getting the fluorescence traces for D1 and D2 MSNs
[dF1,traces1,rawF1,roiList1] = extractROIsBin(myObj,roiGroup,sliceNum,channelNum);
data.dF1.dF = dF1; %DF/F for D1 Neurons
data.dF1.roiList = roiList1; %ROI info for D1 Neurons

roiGroup = 2; %Getting D2 Neuron Data

% Getting the fluorescence traces for D1 and D2 MSNs
[dF2,traces2,rawF2,roiList2] = extractROIsBin(myObj,roiGroup,sliceNum,channelNum);
data.dF2.dF = dF2; %DF/F for D2 Neurons
data.dF2.roiList = roiList2; %ROI info for D2 Neurons

roiGroup = 3; %Getting Other Neuron Data

% Getting the fluorescence traces for ChIs
[dF3,traces3,rawF3,roiList3] = extractROIsBin(myObj,roiGroup,sliceNum,channelNum);
data.dF3.dF = dF3; %DF/F for Other Neurons
data.dF3.roiList = roiList3;%ROI info for Other Neurons

clearvars -except data

%% Behavior & Mirror: Opening and Synchronization
totalTime = 600; % in seconds
timeFreq = 2000; % in hz

firstFile = 4; % Aquitition number CHANGE
lastFile = firstFile;  %%firstFile originale
mirror = true; % This is for the getBehAndMir() function and determines which scanimage file to read
wheelRadius = 0.06; % in m
circumference = 2*wheelRadius*pi; %in m

warning off; % To turn off warnings to speed up code

% Opening behavior and mirror files according to the filenames
[data.rawData, mirArray] = getBehAndMir(firstFile, lastFile, mirror);

% Unwrap the periodic behavior signal to get total distance
data.finalData = unwrapBeh(data.rawData);
data.finalData = circumference*data.finalData;

% Determining onset/offset indices for mirror
mirThreshold = 4; % Voltage threshold for considering "up state" of mirror
minBelowTimeMir = 1; % Mirror must be in down state for at least 1 data point
minAboveTimeMir = 0; % Mirror must be in up state for at least 0 data points
behavior = false; % This is for the getOnsetOffset() function and only impacts how the code handles when the acquisition ends in the middle of a bout or frame

[frameStarts, frameEnds] = ...
    getOnsetOffset(mirArray, mirThreshold, minBelowTimeMir, minAboveTimeMir, behavior);

% Averaging behavior data during each corresponding frame
data.behFrame = [];
for frameInd = 1:length(frameStarts)
    data.behFrame(frameInd) = mean(data.finalData(frameStarts(frameInd):frameEnds(frameInd)));
end

% In case there is one more frame than behavior, cut off final frame
if size(data.dF1.dF, 2) > length(data.behFrame);
    data.dF1.dF = data.dF1.dF(:,1:length(data.behFrame));
    data.dF2.dF = data.dF2.dF(:,1:length(data.behFrame));
    data.dF3.dF = data.dF3.dF(:,1:length(data.behFrame));
else if size(data.dF1.dF, 2) < length(data.behFrame);
        data.behFrame = data.behFrame(1:size(data.dF1.dF, 2));
    end
end

framerate = length(data.behFrame)/totalTime; %in hz
data.time = 1/timeFreq : 1/timeFreq : length(data.finalData)/timeFreq; %time scale
vel = framerate*[diff(data.behFrame), data.behFrame(end) - data.behFrame(end-1)]; %in m/s
data.vel = medfilt1(vel,10); %in m/s

accel = framerate*[diff(data.vel), data.vel(end) - data.vel(end-1)];
data.accel = medfilt1(accel,5);

data.framerate = framerate; %in hz
data.frameTime = (1:length(data.behFrame))/data.framerate;
data.totalDistance = sum(abs(data.vel))/data.framerate;
data.acqNum = firstFile; %Acquisition number for file naming
data.totalTime = totalTime; %Total time elapsed in s
data.timeFreq = timeFreq; %Acquisition rate of behavior in Hz

clearvars -except data
%% Creating Movement Onset/Offset Data

for condition = 1:3; % Run through all conditions (cell types)

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
    behavior = true;
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
    onsetsBeh = iterToMin(signal,onsetsFinal,threshold,true); % Final onsets array (***I don't know if this was included in the paper...)
    offsetsBeh = iterToMin(signal,offsetsFinal,threshold,false); % Final offsets array (***I don't know if this was included in the paper...)

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

% This method of detecting peaks first smooths the calcium trace, then
% finds peaks that are at least 4 frames wide at half prominence, and that
% have a minimum height and prominence of at least 5 times the standard
% deviation of the noise of the calcium trace

% Frequency is calculated by counting the number of peaks within movement
% bouts and during periods of rest, and dividing it by the total amount of
% time spent in those conditions, respectively.

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

% This method of detecting peaks first smooths the calcium trace, then
% finds peaks that are at least 4 frames wide at half prominence, and that
% have a minimum height and prominence of at least 5 times the standard
% deviation of the noise of the calcium trace

% Peak amplitude is calculated by taking the calcium signal at each point
% detected to be a peak

% Rise Time is calculated by determining the time it takes to rise from 20%
% of the peak height to 80% of the peak height

% Decay time is calculated by deterimining the time it takes for the signal
% to decay to 1/3 of the peak height

% Peak width is determined at half of the peak's prominence

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

% This calculates the average fluorescence while the mouse is above or
% below a specific velocity threshold

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

% This saves files according the the mouse name, acquisition date and
% acquisition number

cd('L:\tritsn01labspace\Marta\2p_Data\Calcium Data\')
save([data.imageFile,'_',data.date,'_Acq',num2str(data.acqNum),'.mat'],'data')
cd ..

%%  Plotting Behavioral Data + Calcium Data

% Plot all behavioral data
figure();
h(1)= subplot(4,1,1);
plot(data.time, data.rawData)
title('raw')
ylim([0 1.1])
h(2) = subplot(4,1,2);
plot(data.frameTime, data.behFrame);
title('distance')
h(3)= subplot(4,1,3);
plot(data.frameTime, data.vel)
title('velocity')
h(4)= subplot(4,1,4);
plot(data.frameTime, data.accel)
title('acceleration')
linkaxes(h,'x')
xlabel('Time (s)')

% Uncomment the one you want:
roi = data.dF1.dF;
% roi = data.dF2.dF;
% roi = data.dF3.dF;

% Plot velocity and calcium data for all cells
figure()
for ii =1:size(roi, 1);  %%% CHANGE if you want one cell at time
    title('fluorescence')
    sh(1)= subplot(2,1,2);
    plot(data.frameTime,data.vel)
    title('vel')
    xlabel('Time (s)')
    sh(2)= subplot(2,1,1);
    plot(data.frameTime,roi(ii,:))
    hold on
   end
hold off

linkaxes(sh,'x')


clearvars -except data

%% Plotting at Onset/Offset

% Uncomment the one you want:
roi = data.dF1;
% roi = data.dF2;
% roi = data.dF3;

%frame information (both behavior and imaging)
framerate = data.framerate;
vel = data.vel;
fRatio = size(roi.dF,2)/length(vel);
timeBefore = 4*framerate; % Time before onset/offset (coeff is in seconds)
timeAfter = 4*framerate; % Time after onset/offset (coeff is in seconds)
% scaleFactor = max(data.vel)/max(max(roi.dF)); % Scaling to put DF on same scale as velocity

% Plotting all DF/F and Velocity Traces for ONSET
figure();
subplot(2,1,1)
title('DF/F Traces at Onset')
ylabel('Fluorescence (dF/F)')
hold on
for trace = 1:length(data.indOnsets)
    plot(data.timeDF, roi.onsetsMatrixDF(trace,:));  
end
hold off
subplot(2,1,2)
title('Velocity Traces at Onset')
ylabel('Velocity(m/s)')
xlabel('Time (s)')
hold on
for trace = 1:length(data.indOnsets)
    plot(data.timeDF, roi.onsetsMatrixBeh(trace,:));  
end
hold off

% Plotting all DF/F and Velocity Traces for OFFSET
figure();
subplot(2,1,1);
title('DF/F Traces at Offset')
ylabel('Fluorescence (dF/F)')
hold on
for trace = 1:length(data.indOnsets)
    plot(data.timeDF, roi.offsetsMatrixDF(trace,:));  
end
hold off
subplot(2,1,2);
title('Velocity Traces at Offset')
ylabel('Velocity(m/s)')
xlabel('Time (s)')
hold on
for trace = 1:length(data.indOnsets)
    plot(data.timeDF, roi.offsetsMatrixBeh(trace,:));  
end
hold off

% Plotting average DF/F and Velocity for ONSET and OFFSET
fig = figure;
left_color = [0 0 0];
right_color = [1 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

subplot(2,1,1);
title('Avg DF/F and Velocity at Onset')
yyaxis right
plot(data.timeDF,mean(roi.onsetsMatrixDF,1),'r');
ylabel('Fluorescence (dF/F)')
yyaxis left
plot(data.timeDF,mean(roi.onsetsMatrixBeh,1),'k');
ylabel('Velocity (m/s)')

subplot(2,1,2);
title('Avg DF/F and Velocity at Offset')
xlabel('Time (s)')
yyaxis right
plot(data.timeDF,mean(roi.offsetsMatrixDF,1),'r');
ylabel('Fluorescence (dF/F)')
yyaxis left
plot(data.timeDF,mean(roi.offsetsMatrixBeh,1),'k');
ylabel('Velocity (m/s)')

% Plotting DF/F between ONSET and OFFSET
figure();
for onToOff = 1:length(roi.onsetToOffsetDF)
    subplot(length(roi.onsetToOffsetDF),1,onToOff)
    hold on
    plot(roi.onsetToOffsetTime{onToOff},roi.onsetToOffsetDF{onToOff});
    stem(0,0.3)
    stem(roi.onsetToOffsetTime{onToOff}(end - round((timeBefore + 1)*fRatio)),0.3);
    hold off
    if onToOff == 1
        title('Onset to Offset')
    end
end
xlabel('Time (s)')

% Plotting Velocity with Onsets and Offsets
figure();
title('Velocity with Onsets and Offsets')
ylabel('Velocity(m/s)')
xlabel('Time (s)')
hold on
plot(data.frameTime, vel);
stem(data.indOnsets/framerate,.05*ones(length(data.indOnsets),1))
stem(data.indOffsets/framerate,.05*ones(length(data.indOffsets),1))
hold off

clearvars -except data

%% Plotting Frequency during rest and motion

% Uncomment the one you want:
% roi = data.dF1; 
roi = data.dF2;

figure();
behaviorOrRest=[mean(roi.peakFreqRestMot(:,2)), mean(roi.peakFreqRestMot(:,1))];
bar(1.5, behaviorOrRest(1), 'FaceColor', 'b'); hold on
bar(2.5, behaviorOrRest(2), 'FaceColor', 'r');
title('Number of Peaks During Movement or Rest');
set(gca, 'XTick', [1.5 2.5], 'XTickLabel', {'Movement','Rest'});
ylabel ('Peaks Per Second');

clearvars -except data

%% Plotting all D1/D2 together (can delete)

plottingAllCellsFinal;

clearvars -except data
