function [finalSignal] = fourierFilt(signal,Fc,Fs,varargin)

% [finalSignal] = fourierFilt(signal,Fc,Fs,varargin)
% 
% Summary: This function filters data using the fourier transform and then
% zeroing out undesired frequencies. NOTE: this may create some ringing
% artifacts, and it may be better to just downsample or create a
% butterworth filter. So, caution before using this.
%
% Inputs:
%
% 'signal' - the data to be filtered.
%
% 'Fc' - the cutoff frequency in Hz. This must be a single number for high 
% pass or low pass filters and an array of 2 values for band pass or band 
% stop filters.
% 
% 'Fs' - the sampling rate in Hz.
% 
% 'varargin' - 'high', 'low', 'pass', or 'stop' to determine the type of
% filter (high pass, low pass, band pass, or band stop, specifically).
%
% Outputs:
%
% 'finalSignal' - the filtered output signal. This will be a complex
% number.
%
% Author: Jeffrey March, 2018

T = 1/Fs;             % Sampling period       
L = size(signal,2);   % Length of signal
t = (0:L-1)*T;        % Time vector
finalSignal = zeros(size(signal)); % Initializing final signal vector/matrix

% High Pass Filter
if strcmp(varargin, 'high') && length(Fc) == 1
    for i = 1:size(signal,1)
        X = signal(i,:);
        Y = fft(X);
        Y(1:Fc*L/Fs) = 0;
        Y(L-Fc*L/Fs:end) = 0;
        finalSignal(i,:) = ifft(Y);
    end
end

% Low Pass Filter
if strcmp(varargin, 'low') && length(Fc) == 1
    for i = 1:size(signal,1)
        X = signal(i,:);
        Y = fft(X);
        Y(Fc*L/Fs:L - Fc*L/Fs) = 0;
        finalSignal(i,:) = ifft(Y);
    end
end

% Band Pass Filter
if strcmp(varargin, 'pass') && length(Fc) == 2
    for i = 1:size(signal,1)
        X = signal(i,:);
        Y = fft(X);
        Y(1:Fc(1)*L/Fs) = 0;
        Y(L-Fc(1)*L/Fs:end) = 0;
        Y(Fc(2)*L/Fs:L - Fc(2)*L/Fs) = 0;
        finalSignal(i,:) = ifft(Y);
    end
end

%  Band Stop Filter
if strcmp(varargin, 'stop') && length(Fc) == 2
    for i = 1:size(signal,1)
        X = signal(i,:);
        Y = fft(X);
        Y(Fc(1)*L/Fs:Fc(2)*L/Fs) = 0;
        Y(L-Fc(2)*L/Fs:L-Fc(1)*L/Fs) = 0;
        finalSignal(i,:) = ifft(Y);
    end
end

end

