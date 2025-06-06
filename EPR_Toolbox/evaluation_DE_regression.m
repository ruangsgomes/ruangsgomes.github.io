function [Y_DE_eval] = evaluation_DE_regression (DE, Output_DE, X_eval)

% Output_DE = 
% 
%      bestpop: [1x30 double]
%        besty: [14924x1 double]
%     best_par: [-118.6214 0.2307 1.0240e-007 22.2626 23.0282 -0.7312]
%          rsq: '0.83'
%         RMSE: '18'
%          CoD: [80x1 double]
%          SEE: [80x1 double]

% Author: Guilherme Gomes 

bestpop = Output_DE.bestpop;
best_par = Output_DE.best_par;

m = DE.m; k = DE.k; 

% Compute model from best population and parameters
mat = (reshape(bestpop, m, k));
%mat = (reshape(bestpop, k, m))';

% create auxiliary variable for exponential operation
R = zeros(numel(X_eval(:,1)),k,m);
RM = zeros(numel(X_eval(:,1)),m);
        
% Built model with Bias
if strcmp(DE.bias, 'Yes')
    % evaluate exponential operation for each parameter m
    for j = 1:m
        for kk = 1:k
            R(:,kk,j) = X_eval(:,kk).^mat(j,kk);
        end
        RM(:,j) = prod(R(:,:,j),2) .* best_par(j+1);
    end

    Y_DE_eval = best_par(1) + sum(RM,2);
   
% Built model with NO Bias    
else 
    % evaluate exponential operation for each parameter m
    for j = 1:m
        for kk = 1:k
            R(:,kk,j) = X_eval(:,kk).^mat(j,kk);
        end
        RM(:,j) = prod(R(:,:,j),2) .* best_par(j);%+1);
    end
    %Y_DE_eval = best_par(1) + sum(RM,2);
    Y_DE_eval = sum(RM,2);
end
