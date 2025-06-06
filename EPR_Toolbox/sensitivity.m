function [X_sensitivity, Y_sensitivity] = sensitivity(DE, Output,X)

% This function reads the EPR mathematical structure (DE and Output)
% and the parameter values (X) to automatically run a
% a sensitivity analysis. 

% This function automatically reads the number of parameters analyzed, 
% creates a matrix for each parameter containing a range of values
% between max and min for the analyzed parameter, whereas the non-analyzed
% parameters are keeps at their mean values for all the data-point.
% The matrizes for each parameter is stored (X_sensitivity)
% and also its corresponding output is also stored (Y_sensitivity)

% Date 12.05.2021


% Number of independent variables for sensitivity
n_x = size(X,2);

% Number of data-points
n_dp = size(X,1);
%n_dp = 150;

% Initialize structure for ranges (max, min), and mean values
X_sensitivity.mean = zeros(n_dp,n_x);
X_sensitivity.maxmin = zeros(n_dp,n_x);

% Build struct with mean, and min/max ranges for each variable
for i=1:n_x
    X_sensitivity.mean(:,i) = linspace(mean(X(:,i)),mean(X(:,i)),n_dp);
    X_sensitivity.maxmin(:,i) = linspace(min(X(:,i)),max(X(:,i)),n_dp);
end

% Create big n_x dimension matrix to store X values for sensitivity
X_sensitivity.X = zeros(n_dp,n_x,n_x);

% Fill the big matrix with ranges and mean values
for i=1:n_x
    % Fill matrix of current variable with mean values
    X_sensitivity.X(:,:,i)=X_sensitivity.mean(:,:);
    % Get variable matrix, and replace its column with range values
    X_sensitivity.X(:,i,i)=X_sensitivity.maxmin(:,i);
end

% Create sensitivity output for each variable
Y_sensitivity = zeros(n_dp,n_x);

% Run equation with input matrix of sensitivity
for i=1:n_x
[Y_sensitivity(:,i)] = evaluation_DE_regression (DE, Output, X_sensitivity.X(:,:,i));
end


