function downsampled = downsampleBins(signal, binSize)

% downsampled = downsampleBins(signal, binSize)
% 
% Summary: This function downsamples data by averaging the number of data
% points defined by 'binsize' into one point. It has the advantage of
% reducing noise, while not applying a rolling blurring (smoothing) filter
% to all the data, possibly distorting the time course of the data. It is
% advised to make the binSize evenly divide into the signal length.
%
% Inputs:
%
% 'signal' - the original signal.
%
% 'binSize' - the number of data points you want to average in the
% downsampling.
%
% Outputs:
%
% 'downsampled' - the downsampled signal.
%
% Author: Jeffrey March, 2018

downBins = 0:binSize:length(signal);
if mod(length(signal),binSize) ~= 0
    display('Warning: the bins do not fit evenly into the signal length');
end

for i = 1:length(downBins) - 1
    downsampled(i) = mean(signal(downBins(i) + 1 : downBins(i+1)));
end

end