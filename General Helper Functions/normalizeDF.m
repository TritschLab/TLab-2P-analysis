function normalizedSignal = normalizeDF(signal)

% normalizedSignal = normalizeDF(signal)
%
% Summary:  This function takes a signal, subtracts the baseline (median
%of the entire data set), and then divides everything by the maximum signal
%of that cell.
%
% Inputs:
% 
% 'signal' - the original, full data trace
% 
% Outputs:
% 
% 'normalizedSignal' - the normalized signal
%
% Author: Jeffrey March, 2018

% Making sure signal (and, as a result, the output) is horizontal
if size(signal,1) > size(signal,2);
    signal = signal';
end

normalizedSignal = zeros(size(signal)); % initializing output

% Performing the Normalization
for nCell = 1:size(signal,1);
    baselineAdjust = median(signal(nCell,:));
    baselinedSignal = signal(nCell,:) - baselineAdjust;
    normalizedSignal(nCell,:) = baselinedSignal/max(baselinedSignal);
end

end

