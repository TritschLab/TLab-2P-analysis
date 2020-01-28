%% Summary:
% 
% This script opens and updates all selected files to contain the correct
% rest onset and offset indices.
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

[trials, pathname] = uigetfile('*.mat','MultiSelect','on');
cd(pathname)

if ~iscell(trials)
    tempTrials = trials;
    trials = cell(1);
    trials{1} = tempTrials;
end

% cd('L:\tritsn01labspace\Marta\2p_Data\Calcium Data\sL7\')
% cd('C:\MATLAB\Calcium Data\')
totalFiles = 0;

for trial = 1:length(trials);
    totalFiles = totalFiles + 1               
    load(trials{trial});  
    
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
    onsetsFinal = onsetsFinal + timeShift; % Shifting onset by 0.5 seconds
    offsetsFinal = offsetsFinal - timeShift; % Shifting offset by 0.5 seconds

    data.indOnsetsRest = onsetsFinal;
    data.indOffsetsRest = offsetsFinal;
                
	save(trials{trial},'data')                

end
