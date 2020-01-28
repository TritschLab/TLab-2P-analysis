function [rawData, mirArray] = getBehAndMir(firstFile, lastFile, mirror)

% [rawData, mirArray] = getBehAndMir(firstFile, lastFile, mirror)
%
% Summary: This function concatenates the rotary encoder and mirror data
% and then normalizes rotary encoder data.
%
% Inputs:
%
% 'firstFile' - the first scanimage acquisition trace to be concatenated.
% This is the x in AD0_x.
%
% 'lastFile' - the last scanimage acquisition trace to be concatenated. This
% is the x in AD0_x. This is equal to firstFile if you want to look at just
% one acquisition.
%
% 'mirror' - this is a boolean value: true if you're acquiring mirror data,
% and false if you're not. The function concatenates these mirror
% acquisitions in the same way as the behavior data.
%
% Outputs:
%
% 'rawData' - the normalized, concatenated rotary encoder data. Values range
% from 0-1 rotation.
%
% 'mirArray' - the concatenated mirror data. [] if mirror == false.
%
% Author: Jeffrey March, 2018

behArray = [];
for fileNameIndex = firstFile:lastFile
    file = getfield(open(sprintf('AD0_%d.mat',fileNameIndex)), sprintf('AD0_%d',fileNameIndex));
    behArray = cat(2,behArray,file.data);
end
rawData = behArray;
rawData = (rawData - min(rawData)); % Making minimum point = 0
rawData = rawData/max(rawData); % Max point = 1

mirArray = [];
if mirror == true
    for fileNameIndex = firstFile:lastFile
        mirfile = getfield(open(sprintf('AD1_%d.mat',fileNameIndex)), sprintf('AD1_%d',fileNameIndex));
        mirArray = cat(2,mirArray,mirfile.data);
    end
end

end

