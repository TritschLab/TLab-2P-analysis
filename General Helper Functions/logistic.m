function y = logistic(k,x)

% y = logistic(k,x)
% 
% Summary: This is a simple logistic function, used for the purposes of
% fitting data.
%
% Inputs:
%
% 'k' - a parameter.
%
% 'x' - an array of x values.
%
% Outputs:
%
% 'y' - the logistic output.
%
% Author: Jeffrey March, 2018

y = zeros(size(x));
for i = 1:length(x)
    y(i) = 2*(-0.5 + 1/(1 + exp(-k*x(i))));
end

end