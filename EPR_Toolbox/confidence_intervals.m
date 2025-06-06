function [low, up] = confidence_intervals (simulation)
% sort FS simulations to get confidence interval
s = sort(simulation);

% number of posterior solutions
np = size(s,1);

for i = 1:size(s,2)
    % get the column
    col = s(:,i);
    % 0.975
    low(:,i) = col(ceil(0.025*np));
    up(:,i) = col(ceil(0.975*np));
end
