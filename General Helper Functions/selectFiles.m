function [trials, path] = selectFiles()

% [trials, path] = selectFiles()
%
% Summary: This function runs the uigetfile() function, changes the path,
% and makes sure, if one trial is selected, that it is placed in a cell
% array.
%
% Inputs:
%
% User-selected .mat file names
% 
% Outputs:
%
% 'trials' - user-selected .mat file names in a cell array.
%
% 'path' -  a string containing the path of the user-selected filenames.
%
% Author: Jeffrey March, 2018

[trials, path] = uigetfile('*.mat','MultiSelect','on');

cd(path)

if ~iscell(trials)
    tempTrials = trials;
    trials = cell(1);
    trials{1} = tempTrials;
end

end