function [Y_cal, X_cal, Y_eval, X_eval] = split_dataset(Y,X,P)

for i=1:20
    [m,n] = size(Y) ; % get lines (m) and column (n) sizes
    %P = 0.75 ;% 75% for calibration
    idx = randperm(m)  ;
    Y_cal_merged(:,i) = Y(idx(1:round(P*m)),:) ; 
    Y_eval_merged(:,i) = Y(idx(round(P*m)+1:end),:) ;
    % Loop to store independent variables into 3D matrix
    for j=1:size(X,2)
         X_cal_merged(:,i,j) = X(idx(1:round(P*m)),j) ; 
         X_eval_merged(:,i,j) = X(idx(round(P*m)+1:end),j) ;
    end
    
    %Compute mean and standard deviation of dependent variables
    Y_cal_mean = mean(Y_cal_merged(:,i));
    Y_eval_mean = mean(Y_eval_merged(:,i));
    Y_cal_std = std(Y_cal_merged(:,i));
    Y_eval_std = std(Y_eval_merged(:,i));
    
    % Compute distance between mean and standard dev
    Dist_mean(i) = (Y_cal_mean - Y_eval_mean)^2;
    Dist_std(i) = (Y_cal_std - Y_eval_std)^2;
end

% Now rank them first for mean
% Rank the models considering all parameters and statistics
[sd_BestModel,ind_BestModel]=sort(Dist_mean,'ascend'); % sort function
Best_Model =ind_BestModel(1); %the corresponding indices

% Compute best model
Y_cal = Y_cal_merged(:,Best_Model);
Y_eval = Y_eval_merged(:,Best_Model);

% Compute best independents
for j=1:size(X,2)
    X_cal(:,j) = X_cal_merged(:,Best_Model,j);
    X_eval(:,j) = X_eval_merged(:,Best_Model,j);
end
