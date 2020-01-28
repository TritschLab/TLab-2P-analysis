function [sliceNum, channelNum] = getDF(myObj)

% [sliceNum, channelNum] = getDF(myObj)
%
% Summary: This function opens the Acuqisition2P object and prepares it for
% extraction of ROI information. Used in martaCleanCodeUIFINAL.
%
% Inputs:
%
% 'myObj' - the Acquisition2P object.
%
% Outputs:
%
% 'sliceNum' - the downsampled signal.
%
% 'channelNum' - the downsampled signal.
%
% Author: Jeffrey March, 2018

movNum = 2;
castType = 'single';
sliceNum = 1;
channelNum = myObj.motionRefChannel;
mov = readCor(myObj,movNum,castType,sliceNum,channelNum);
rawMov = readRaw(myObj,movNum,castType);

sliceNum = 1; %Choose a slice to analyze
channelNum = 1; %Choose the GCaMP channel
movNums = []; %this will default to all movie files
radiusPxCov = 11; %default, may need zoom level adjustment
temporalBin = 8; %default (tested w/ 15-30hz imaging), may need adjustment based on frame rate
writeDir = []; %empty defaults to the directory the object is saved in (the 'defaultDir')

img = myObj.meanRef;
img(img<0) = 0;
img(isnan(img)) = 0;
img = sqrt(img);
img = adapthisteq(img/max(img(:)));

actImg = myObj.roiInfo.slice(sliceNum).covFile.activityImg;
% img = img/2 + actImg/2;

smoothWindow = 15; % Gaussian window with std = smoothWin/5, for displaying traces
excludeFrames = []; %List of frames that need to be excluded, e.g. if they contain artifacts
myObj.selectROIs(img,sliceNum,channelNum,smoothWindow,excludeFrames);

% Take the time to read through 'Using the ROI selection tool', and then
% select your cells. Once you've selected ROIs, be sure to save the acquisition again
myObj.save;

end

