function plotAllRows(data)

% plotAllRows(data)
% 
% Summary: This function takes a 2D matrix where each row is a data trace 
% and plots all rows superimposed on one another.
%
% Inputs:
%
% 'data' - a 2D matrix, where each row is a data trace.
%
% Outputs:
%
% a figure containing all rows plotted, superimposed on one another.
%
% Author: Jeffrey March, 2018

figure;
hold on
for i = 1:size(data,1)
    plot(data(i,:));
end

end