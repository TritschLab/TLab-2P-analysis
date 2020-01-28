%% Summary:
% 
% This script is for running multiple image registrations overnight for the
% Acquisition2P objects to be available for manual cell selection. The last
% section of the script can be used to move the objects from one directory
% to another. Instructions are given below.
% 
% Inputs:
% 
% uncorrected movie files
% 
% 'myJobDir' - user-defined directory containing unprocessed acq2p objects
% and uncorrected movie files.
% 
% 'destDir' - destination folder for transferring objects into a new
% directory.
%
% Outputs:
% 
% Acquisition2P object
% 
% Author: Jeffrey March, 2018

%% Initialization

% First initialize objects into each Folder. The code will prompt you to
% select files, and you will highlight all the movie files in one folder,
% thus leading to the creation of the Acquisition2P object.
% NOTE: You must copy all the uncorrected movies into that folder as well.

myObj = Acquisition2P([],@SC2Pinit);
myObj.save

%% Doing the job for each folder (Run Overnight)

% Can have multiple objects per date with DIFFERENT object names... ie.
% SW024 acq 1-2, SW025 acq 1-2, and SW026 acq 1-2 can all be in the same
% folder. Then run acq2pJob Processor on that folder. Or they can all be in
% different folders. You'll just need to copy and paste the three line
% piece of code again and again...

myJobDir = 'L:\tritsn01labspace\Marta\2p_Data\SW064\180130 acq3\'; % User-defined directory containing unprocessed acq2p objects.
ajp = acq2pJobProcessor(myJobDir);

clear all

%% Changing directories

% Use the newDir function to move objects into a new directory, defined by
% destDir.

destDir = 'C:\MATLAB\Marta\SW064\180131 acq4';
myObj = SW064;
newDir(myObj,destDir,0,1,1);