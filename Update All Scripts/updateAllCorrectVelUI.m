%% Summary:
% 
% This script opens and updates all selected files to contain the correct
% velocity information.
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

[trials] = uigetfile('*.mat','MultiSelect','on');

% cd('C:\MATLAB\Calcium Data\')
totalFiles = 0;
scale = 1.5;

for trial = 1:length(trials);
	load(trials{trial});
	totalFiles = totalFiles + 1
    
    data.finalData = scale*data.finalData;
    data.behFrame = scale*data.behFrame;
    data.vel = scale*data.vel;
    data.accel = scale*data.accel;
    data.totalDistance = scale*data.totalDistance;
    
    data.dF1.onsetsMatrixBeh = scale*data.dF1.onsetsMatrixBeh;
    data.dF1.offsetsMatrixBeh = scale*data.dF1.offsetsMatrixBeh;
    data.dF1.meanMovVel = scale*data.dF1.meanMovVel;
    data.dF1.meanRestVel = scale*data.dF1.meanRestVel;
    data.dF1.meanTotVel = scale*data.dF1.meanTotVel;
    data.dF1.meanBinVel = scale*data.dF1.meanBinVel;

    
    data.dF2.onsetsMatrixBeh = scale*data.dF2.onsetsMatrixBeh;
    data.dF2.offsetsMatrixBeh = scale*data.dF2.offsetsMatrixBeh;
    data.dF2.meanMovVel = scale*data.dF2.meanMovVel;
    data.dF2.meanRestVel = scale*data.dF2.meanRestVel;
    data.dF2.meanTotVel = scale*data.dF2.meanTotVel;
    data.dF2.meanBinVel = scale*data.dF2.meanBinVel;
    
    for i = 1:length(data.dF1.onsetToOffsetBeh)
        data.dF1.onsetToOffsetBeh{i} = scale*data.dF1.onsetToOffsetBeh{i};
        data.dF2.onsetToOffsetBeh{i} = scale*data.dF2.onsetToOffsetBeh{i};
    end
	save(trials{trial},'data')
end
cd ..