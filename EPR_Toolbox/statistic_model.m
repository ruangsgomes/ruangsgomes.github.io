function [Out_stat] = statistic_model (obs, model) 

% R-squared statistic (evaluation)
residual = obs - model;
SSresid = sum(residual.^2);
SStotal = (length(obs-1) * var(obs));
rsq = 1 - SSresid/SStotal;
%rsq = num2str(rsq,2);

% RMSE
RMSE = sqrt(SSresid/numel(obs));
%RMSE = num2str(RMSE,2);

SSE = SSresid/numel(obs);

% coef determinacao Shahin ou Nash and Sutcliffe
den = sum((obs - mean(obs)).^2);
R2 = 1 - SSresid/den;

% MAE Shahin
MAE = sum(abs(residual))/numel(obs);

% r Shahin
r = sum((obs - mean(obs)).*(model - mean(model)))/...
    sqrt(sum((obs - mean(obs)).^2)*sum((model - mean(model)).^2));

% Phogat et al 2016
% Statistical Assessment of a Numerical Model Simulating 
% Agro Hydrochemical Processes in Soil under Drip Fertigated 
% Mandarin Tree
% Relative efficiency
E_rel = 1 - (sum(((obs - model)./(obs)).^2)/...
    sum(((obs - mean(obs))./(mean(obs))).^2));

% Percentage Bias
PBIAS = 100*sum(model - obs)/sum(obs);

%% output
Out_stat.rsq = rsq;
Out_stat.RMSE = RMSE;
Out_stat.SSE = SSE;
Out_stat.R2 = R2;
Out_stat.MAE = MAE;
Out_stat.r = r;
Out_stat.E_rel = E_rel;
Out_stat.PBIAS = PBIAS; 