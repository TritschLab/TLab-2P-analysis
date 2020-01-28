%% Summary:
% 
% This script opens and updates all selected files to contain the correct
% total distance traveled for a given acquisition.
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

cd('C:\Users\labadmin\Documents\MATLAB\Calcium Data\')
totalFiles = 0;
for trial = 1:length(trials);
	load(trials{trial});
	totalFiles = totalFiles + 1
	data.totalDistance = sum(abs(data.vel))/data.framerate;
%	 data.totalDistance = sum(abs(data.vel))/data.sampleRate;
	data = rmfield(data,'totaldistance');

	save(trials{trial},'data')
end
cd ..