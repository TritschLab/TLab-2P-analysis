function [finalData, flagMatrix] = unwrapBeh(rawData)

% [finalData, flagMatrix] = unwrapBeh(rawData)
%
% Summary:  This function converts the periodic encoder trace to a smooth
% measure of number of wheel rotations. It can work for any periodic signal
% that goes from 0 to 1. Should be robust to extra data points found *at*
% the period threshold, in between the max and min.
%
% Inputs:
%
% 'rawData' - the periodic, rotary encoder trace. This trace should be
% normalized, having a min value of 0 and max value of 1.
%
% Outputs:
%
% 'finalData' - the final unwrapped data that will be a smooth measure of
% number of periods (usually, wheel rotations) in the signal.
%
% 'flagMatrix' - the corresponding flag for each data point. This is mostly
% for debugging purposes.
%
% Author: Jeffrey March, 2018

% Initializing variables for loop
lapCounter = 0; % This keeps track of number of laps
index = 2; % Initializing index for 'while' loop
finalData = zeros(1,length(rawData)); % Initializing finalData
lapCrossings = zeros(1,length(rawData)); % Initializing lapCrossings
flagMatrix = zeros(1,length(rawData)); % Initializing flagMatrix
finalData(1) = rawData(1); % Defining index = 1 point for finalData
lapCrossings(1) = 0; % Defining index = 1 point for lapCrossings

% Creating initial flag values
if rawData(1) >= 0.45
    flag = 1;
else if rawData(1) < 0.45
        flag = 0;
    end
end
flagMatrix(1) = flag;

% Looping through data and connecting successive laps
while index <= length(rawData)
    
    if rawData(index - 1) - rawData(index)  > 0.4 && flag == 1 % Checking if a lap is complete
        rawData(index) = floor(rawData(index)); % After crossing the the period threshold forward, the next point gets the value 0
        lapCounter = lapCounter + 1; % Adding one lap
        lapCrossings(index) = 1; % Value of 1, signifying lap crossing
        flag = 0;
    elseif rawData(index) - rawData(index - 1) > 0.4 && flag == 0 % Checking if mouse has gone backwards at lap crossing
        rawData(index) = ceil(rawData(index)); % If going backward across the period threshold, the next point gets the value 1
        lapCounter = lapCounter - 1; % Subtracting one lap
        lapCrossings(index) = -1; % Value of -1, signifying backwards lap crossing
        flag = 1;
        
    % Account for if mouse runs backwards at lapCounter threshold (This was
    % previously before the other if/elseif statements)
    elseif rawData(index) >= 0.45 && rawData(index) < 0.55
        flag = 1;
    elseif rawData(index) >= 0.35 && rawData(index) < 0.45
        flag = 0;
    end
    
    finalData(index) = rawData(index) + lapCounter; % Updating the index by adding the total number of laps
    flagMatrix(index) = flag;
    
    if lapCrossings(index) == 0;
        lapCrossings(index) = 0; % Value of 0, signifying no lap crossed
    end
    
    index = index + 1;
end

end

