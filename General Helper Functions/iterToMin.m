function [localMinsFinal] = iterToMin(signal, array, threshold, isOnset)

% [localMinsFinal] = iterToMin(signal, array, threshold, isOnset)
%
% Summary: This function iterates towards local minimum values,
% sequentially from a given array of indices
%
% Inputs:
%
% 'signal' - the signal in which we are looking for minimums
% 
% 'array' - the indices of the starting points
%
% 'threshold' - the minimum to iterate toward
%
% 'isOnset' - if true, the iteration goes left, if false, the iteration
% goes right (this is a holdover from using this code for onsets and
% offsets of mouse movement bouts)
% 
% Outputs:
%
% 'localMinsFinal' - an array of the new local minimums, iterated towards
% from the initial 'array' 
%
% Author: Jeffrey March, 2018

localMinsFinal = zeros(size(array)); % initialiaing results array

for i = 1:length(array)
    localMin = array(i);
    threshInd = i;
    
    if length(threshold) == 1
        threshInd = 1;
    end
    
    % Iterating towards the local minimum (direction depends on isOnset)   
    while signal(localMin) > threshold(threshInd)
        localMin = localMin - (isOnset*2 - 1);
        
        % Checking to make sure onset/offset doesn't run off end of signal
        if localMin < 1 || localMin > length(signal)
            localMin = nan;
            break
        end
        
    end
    
    localMinsFinal(i) = localMin;
    
end

end

