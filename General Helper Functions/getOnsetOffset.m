function [onsetsFinal, offsetsFinal] = getOnsetOffset(signal, ampThreshold, minBelowTime, minAboveTime, behFlag)

% [onsetsFinal, offsetsFinal] = getOnsetOffset(signal, ampThreshold, minBelowTime, minAboveTime, behFlag)
%
% Summary:  This function determines the indices of onsets and offsets
% for a movement bout or a periodic mirror position.
%
% Inputs:
% 
% 'signal' - the original, full data trace
%
% 'ampThreshold' - the minimum height of the signal to signify an onset
%
% 'minBelowTime' - the minimum number of points that the signal must be
% below the threshold to consider a bout over
%
% 'minAboveTime' - the minimum number of points that the signal must be
% above the threshold to be considered a bout
%
% 'behFlag' - a boolean that determines whether behavioral bouts are being
% computed. A value of 'true' means that if the acquisition ends during a
% bout, the final onset should be eliminated. A value of 'false' means that
% even though the acquisition ends during a bout, the final point will be
% considered an offset.
% 
% Outputs:
% 
% 'onsetsFinal' - the indices of the beginnings of periods of signal
% increase.
% 
% 'offsetsFinal' - the indices of the ends of periods of signal increase.
%
% Author: Jeffrey March, 2018 

aboveThresh = find(signal >= ampThreshold);
onsets = cat(2,aboveThresh(1),aboveThresh(find(diff(aboveThresh) > minBelowTime) + 1));
offsets = cat(2,aboveThresh(find(diff(aboveThresh) > minBelowTime)),aboveThresh(end));
if behFlag == true
    onsets = onsets(1:length(offsets));
end
if behFlag == false
    if length(offsets) < length(onsets)
        offsets = cat(2,offsets,aboveThresh(end));
    end
end
offsetsFinal = offsets((offsets - onsets) > minAboveTime);
onsetsFinal = onsets((offsets - onsets) > minAboveTime);

end

