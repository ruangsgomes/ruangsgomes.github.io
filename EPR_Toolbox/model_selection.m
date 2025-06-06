function [BEST]=model_selection(n,X,Y, X_cal, Y_cal, names)

% This function selects the model that best fits the mean sensitivities
% The selection considers the best model considering all indep. variables
% Which is ranked based on the sensitivity analyses

% Finally all sensitivity graphs are plotted

% 30.05.2021

% Fontsize
fax = 18;

% compute independent variables
n_x = size(X.mean,2);

%% PREPARE DATA

% COMPUTE MEAN VALUES
% Number of parameters
for i=1:n_x
    % Number of MonteCarlo runs
    for j = 1: n 
        Y_mean(i).par(:,j) = Y(j).Y_sensitivity(:,i);
    end
   % Mean values of parameter i
   mean_sensitivity(:,i)=mean(Y_mean(i).par(:,:),2); 
end

% COMPUTE STATISTICS - DISTANCE FROM EACH MODEL AGAINST THE MEAN
% Number of parameters
for i=1:n_x
    
    % Number of MonteCarlo runs
    for j = 1: n 

        % Compute statistics for model j
        [Out_stats] = statistic_model (mean_sensitivity(:,i), Y(j).Y_sensitivity(:,i)) ;
        
        % Check negative values (only rsq, r and E_rel) and replace 
        if Out_stats.rsq < 0
            Out_stats.rsq = 0.01;
        end
        
        if Out_stats.r < 0
            Out_stats.r = 0.01;
        end
        
        if Out_stats.E_rel < 0
            Out_stats.E_rel = 0.01;
        end
        
        % Save statistic for all monte carlo models (j) and ind. variable
        model_stats(i).rsq(j) = Out_stats.rsq;
        model_stats(i).RMSE(j) = Out_stats.RMSE;
        model_stats(i).r(j) = Out_stats.r;
        model_stats(i).E_rel(j) = Out_stats.E_rel;    
    end
end

% COMPUTE SUMMARY STATISTICS
% Number of parameters
for i=1:n_x
    model_stats(i).max.rsq = max(model_stats(i).rsq);
    model_stats(i).max.RMSE = max(model_stats(i).RMSE);
    
    model_stats(i).min.rsq = min(model_stats(i).rsq);
    model_stats(i).min.RMSE = min(model_stats(i).RMSE);

    model_stats(i).mean.rsq = mean(model_stats(i).rsq);
    model_stats(i).mean.RMSE = mean(model_stats(i).RMSE);   
    
    model_stats(i).stdDev.rsq = std(model_stats(i).rsq);
    model_stats(i).stdDev.RMSE = std(model_stats(i).RMSE);
end

%% MODEL SELECTION

% auxiliary variables to compute 
aux_sd_rsq = zeros(n,1,n_x);
aux_ind_rsq = zeros(n,1,n_x);

% initiate variable that contains the indices of the best models
model_rsq = zeros(n,1);
model_RMSE = zeros(n,1);

for i=1:n_x
    
    % Sort descending for R2
    model_sort = model_stats(i).rsq;% Another auxiliary for sort function
    [sd_rsq,ind_rsq]=sort(model_sort,'descend'); % sort function
    aux_sd_rsq(:,:,i)=sd_rsq';  %sorted data
    aux_ind_rsq(:,:,i) =ind_rsq'; %the corresponding indices
    model_rsq = [model_rsq , aux_ind_rsq(:,:,i)];% Merge the ranking of each variable
   
    % Sort descending for RMSE
    model_sort = model_stats(i).RMSE;% Another auxiliary for sort function
    [sd_RMSE,ind_RMSE]=sort(model_sort,'ascend'); % sort function
    aux_sd_RMSE(:,:,i)=sd_RMSE';  %sorted data
    aux_ind_RMSE(:,:,i) =ind_RMSE'; %the corresponding indices
    model_RMSE = [model_RMSE , aux_ind_RMSE(:,:,i)];% Merge the ranking of each variable
       
end

% Exclude the first column (zero values)
model_rsq(:,1) = []; 
model_RMSE(:,1) = []; 


% EVALUATE THE MOST CONSISTENT MODEL OF ALL INPUT PARAMETERS
% AND RETURN ITS INDICE

[BEST] = model_consistent(n,model_rsq, model_RMSE);

%% TEXT DATA

 % output file of the best simulation
    fileID = fopen('Model_Selection.txt','a');
    fprintf(fileID,'%12s\r\n','---------- Model selection output file ----------');
    fprintf(fileID,'%12s\r\n','algorithm: MonteCarlo MODEGA with optimum length');
    fprintf(fileID,'%12s\r\n','Summary statistics of sensitivity');
    for i =1:n_x
    fprintf(fileID,'%12s\r\n','--------------------------------------------');
    fprintf(fileID,'%12s %6.4f\r\n','Parameter =',i);
    fprintf(fileID,'%12s\r\n','--------- R^2 ---------');
    fprintf(fileID,'%12s %6.4f\r\n','Max=', model_stats(i).max.rsq);
    fprintf(fileID,'%12s %6.4f\r\n','Min=', model_stats(i).min.rsq);
    fprintf(fileID,'%12s %6.4f\r\n','Mean=', model_stats(i).mean.rsq);
    fprintf(fileID,'%12s %6.4f\r\n','Std. dev.=', model_stats(i).stdDev.rsq);
    fprintf(fileID,'%12s\r\n','--------- RMSE ---------');
    fprintf(fileID,'%12s %6.4f\r\n','Max=', model_stats(i).max.RMSE);
    fprintf(fileID,'%12s %6.4f\r\n','Min=', model_stats(i).min.RMSE);
    fprintf(fileID,'%12s %6.4f\r\n','Mean=', model_stats(i).mean.RMSE);
    fprintf(fileID,'%12s %6.4f\r\n','Std. dev.=', model_stats(i).stdDev.RMSE);
    fprintf(fileID,'%12s\r\n','--------- Best Model ---------');
    fprintf(fileID,'%12s %6.4f\r\n','R^2=', model_stats(i).rsq(1,BEST));
    fprintf(fileID,'%12s %6.4f\r\n','RMSE=', model_stats(i).RMSE(1,BEST));
    end
    fclose(fileID);

end


function [Best_Model] = model_consistent(n,model_rsq, model_RMSE)
% This model evaluates the most consistent model considering
% all the input parameters and merges the ranking with R2 and RMSE
% statistical indicators

% Number of independente variables (X)
n_x = size(model_rsq,2);

% Number of monte carlo (n)

% Count rsq
count_rsq = zeros(n,n_x);
count_RMSE = zeros(n,n_x);

% THIS LOOP COMPUTES THE POSITION OF EACH INDICE
for i = 1:n_x
    for j = 1:n
       % Get indice number
       aux_rsq = model_rsq(j,i);
       aux_RMSE = model_RMSE(j,i);
       % Sum its position for ranking
       count_rsq(aux_rsq,i) = count_rsq(aux_rsq,i)+j;
       count_RMSE(aux_RMSE,i) = count_RMSE(aux_RMSE,i)+j;
       
    end
end

% Sum all parameters values
count_all_rsq = sum(count_rsq,2);
count_all_RMSE = sum(count_RMSE,2);

% Merge ranking 
Merge_model_aux = [count_all_rsq count_all_RMSE];
Merge_model = sum(Merge_model_aux,2);

% Rank the models considering all parameters and statistics
[sd_BestModel,ind_BestModel]=sort(Merge_model,'ascend'); % sort function
Best_Model =ind_BestModel(1); %the corresponding indices
%Middle_Model = ind_BestModel(n*0.93);
%Worst_Model = ind_BestModel(n*0.98);

end


