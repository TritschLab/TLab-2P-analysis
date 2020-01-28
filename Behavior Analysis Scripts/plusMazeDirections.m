%% Summary:
% 
% This script describes and executes a procedure for automatedly analyzing
% plus maze videos for determining a mouse's turning bias.
% 
% Instructions:
% 
% Draw ROIs in imageJ for each arm, up until the black tape, plus one
% triangle of an area off the pluz maze for background signal variation
% subtraction. Once you get the traces of each arm's average pixel
% intensity over time, copy and paste all five columns into a variable
% called unnamed. You may want to delete the first ten frames or so because
% of weird camera behavior.
% 
% The code basically works by calculating the mean pixel intensity in each
% arm, and when the mouse enters the arm, since the mouse is black, the
% mean pixel intensity of that arm will drop. The code then normalizes and 
% thresholds the mean pixel intensity for each arm at a certain threshold
% to determine when the mouse is in each arm. The code then uses that
% information to determine the order of arm entries, and thus the direction
% of each turn.
% 
% Bugs may arise if the movies have automatic contrasting, which may lead
% to funky shifts in mean pixel intensity. To avoid this, look at some
% movies/traces manually, and if necessary, delete some of the first frames
% from the movie. Another debugging option is to play around with the
% threshold.
% 
% Properties for paper: 
% 
% threshold = 0.75%
% 
% Inputs:
% 
% 'unnamed' - a matrix of mean pixel intensity traces (columns), measured
% in imageJ and copied and pasted into this variable
%
% Outputs:
% 
% 'leftTurns' - the number of left turns the mouse makes
% 
% 'rightTurns' - the number of right turns the mouse makes
% 
% 'straightTurns' - the number of straight turns the mouse makes
% 
% 'directionsWords' - the order of turns in a given acquisition
% 
% Author: Jeffrey March, 2018

%% Run this first

SWplus = unnamed';

%% Normalizing Arms

for i = 1:4
    SWplus(i,:) = SWplus(i,:) - SWplus(5,:);
%     abs((max(SWplus(i,:)) - min(SWplus(i,:)))/max(SWplus(i,:)))
    if max(SWplus(i,:)) - min(SWplus(i,:)) > 2 % Checking to make sure mouse enters each arm; may need to change the hardcoded value
        SWplus(i,:) = SWplus(i,:) - min(SWplus(i,:));
        SWplus(i,:) = SWplus(i,:)/max(SWplus(i,50:end));
    else
        warning('One arm was not crossed! Check manually.')
        SWplus(i,:) = ones(size(SWplus(i,:)));
    end
end

%% Thresholding SW032

threshold = 0.75;
for i = 1:4
    SWplus(i,:) = SWplus(i,:) < threshold;
end

%% Determining which arm the mouse is in at any given time

summed = SWplus(1,:) + 10*SWplus(2,:) + 100*SWplus(3,:) + 1000*SWplus(4,:);

for index = 1:length(summed)
    if mod(log10(summed(index)),1) ~= 0 % Is value == 1,10,100,1000 exactly (not more than one arm at once)
        summed(index) = 0;
    end
    if summed(index) > 0
        summed(index) = log10(summed(index)) + 1; % Making arms 1,2,3,4 instead of 0,1,2,3
    end
end

%% Taking the arms the mouse is in and determining direction changes

armSequence = [];
counter = 1;
for index = 1:length(summed)
    if summed(index) ~= 0
        currentArm = summed(index);
        if counter == 1
            armSequence(counter) = currentArm;
            counter = counter + 1;
        else if currentArm ~= armSequence(counter - 1);
                armSequence(counter) = currentArm;
                counter = counter + 1;
            end
        end       
    end
end

directions = diff(armSequence);
%% Counting left, right, straight turns

% Key: 1 = right; 2 = straight; 3 = left; -1 = left; -2 = straight; -3 = right

leftTurns = length(directions(directions == 3 | directions == -1));
rightTurns = length(directions(directions == 1 | directions == -3));
straightTurns = length(directions(directions == 2 | directions == -2));

%% Translating directions into strings
directionsWords = {};
for i = 1:length(directions)
    switch directions(i)
        case {-1,3}
            directionsWords{i} = 'left';
        case {1, -3}
            directionsWords{i} = 'right';
        case {2, -2}
            directionsWords{i} = 'straight';
    end
end
