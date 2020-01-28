%% Summary:
% 
% This script opens and updates all selected files to contain the correct
% rest onset and offset indices for motor acquisitions, taking into acount
% the video subtraction method of detecting motion.
% 
% Inputs:
% 
% User-selected .mat file names
%
% Outputs:
% 
% Updated data structure, re-saved under its original name in the .mat file
% 
% Author: Jeffrey March, 2018

%% Main Code

clear all

[trials, pathname] = uigetfile('*.mat','MultiSelect','on');

if ~iscell(trials)
    tempTrials = trials;
    trials = cell(1);
    trials{1} = tempTrials;
end

cd(pathname)
totalFiles = 0;

for trial = 1:length(trials);
	load(trials{trial});                               
	totalFiles = totalFiles + 1
    if isfield(data, 'subVideo')
        
        %% Normalizing subVideo
        
        sortedSignal = sort(data.subVideo); % sort signal to get percentile data
        ninety = sortedSignal(round(.9*length(sortedSignal))); % 90th percentile of signal
        ten = sortedSignal(round(.1*length(sortedSignal))); % 10th percentile of signal
        subVideoNorm = (data.subVideo - ten)/(ninety - ten); % 90th percentile = 1, 10th percentile = 0
        
        subVideoNorm = smooth(subVideoNorm,5)';
        
        % figure;
        % plot(subVideoNorm);
        
        %% Getting rest onset and offset for subVideo
        
        framerate = data.framerate;
        timeThreshold = 4*framerate; %CHANGE number 4,5,6,7.. a seconda dimensione finestra pre post onset
        velThreshold = 0.2; %CHANGE was .004
        minRestTime = 4*framerate; %CHANGE number 4,5,6,7.. a seconda dimensione finestra pre post onset
        minRunTime =   1; %CHANGE
        behavior = true; %CHANGE
        signal = subVideoNorm; %This can be either absolute value or not
        timeShift = round(1*framerate);
        
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
        
        data.subVideoRestOnsets = onsetsFinal;
        data.subVideoRestOffsets = offsetsFinal;
        
        %% Plotting Example Trace With Rest Onsets
        
        figure;
        plot(subVideoNorm);
        hold on
        stem(data.subVideoRestOnsets,subVideoNorm(data.subVideoRestOnsets));
        stem(data.subVideoRestOffsets,subVideoNorm(data.subVideoRestOffsets));
        
        %% Saving Data
        
        save(trials{trial},'data')
    end

end

