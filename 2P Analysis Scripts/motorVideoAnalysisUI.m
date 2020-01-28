%% Summary:
% 
% This script calculates and saves the video pixel subtractions for motor
% trials, in order to ultimately calculate "true rest." The process for
% this entails creating a folder called "Videos" in the root directory
% where all the folders of the data structures for each experimental
% condition are located. (E.g. C:\MATLAB\Calcium Data\Lesion FINAL\Videos\mPre
% with C:\MATLAB\Calcium Data\Lesion FINAL\mPre\ etc.). The user selects
% the video files, and the code performs the analysis, taking each frame
% and subtracting it from the previous frame and adding up the absolute
% values of the pixel differences. This creates a trace measuring movement
% aside from that detected by the rotary encoder.
% 
% Inputs:
% 
% **NOTE THE DIFFERENCE** User-selected .avi file names
%
% Outputs:
% 
% 'data.subVideo' - a matrix containing the video pixel subtraction trace,
% updated and saved to the original data structure
% 
% Author: Jeffrey March, 2018

%%
clear all
[trials, path] = uigetfile('*.avi','MultiSelect','on');

if ~iscell(trials)
    tempTrials = trials;
    trials = cell(1);
    trials{1} = tempTrials;
end

totalFiles = 0;
for trial = 1:length(trials);
    tic
    totalFiles = totalFiles + 1
    cd(path)
	file = trials{trial}; % storing the file name for loading purposes
    srcvid = VideoReader(file); % loading the video file                           
	

    L = srcvid.NumberOfFrames; % automatically-detected number of frames
    H = srcvid.Height; % automatically-detected height of frames
    W = srcvid.Width; % automatically-detected width of frames
    diffFrame = zeros(1,L-1);
    resultFrames = zeros(H,W,2);
    for frame = 1:L
        im = read(srcvid,frame);
        if frame == 1
            resultFrames(:,:,1) = rgb2gray(im); % setting the first frame as the previous frame
            resultFrames(:,:,2) = rgb2gray(im); % setting the first frames as the current frame
        else
            resultFrames(:,:,1) = resultFrames(:,:,2); % setting the now previous frame as such
            resultFrames(:,:,2) = rgb2gray(im); % loading the current frame and setting it as such
            diffFrame(frame - 1) = sum(sum(abs(double(resultFrames(:,:,1)) - double(resultFrames(:,:,2))))); % Subtracting successive frames
        end
    end
    
    cd([path(1:36),path(44:end)]); % changing to path of data structure
    load([file(1:end-4),'.mat']) % loading the original data structure
    data.subVideo = diffFrame; % updating the original data structure
    save([file(1:end-4),'.mat'],'data') % saving the updated data structure
    toc
end