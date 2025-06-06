function [Output_MODE] = MODE_algorithm (DE, Y_cal, X_cal)

%% MODE - Multi-objective Differential evolution algorithm

% dimension of the problem
DE.d = DE.m*DE.k;  
% R-matrix for index     
for ii = 1:DE.N, DE.R(ii,1:DE.N-1) = setdiff(1:DE.N,ii); end 
% Parameter information
Par_info.min = zeros(DE.d,1)'; Par_info.min(1:end) = DE.exponents(1);  
Par_info.max = zeros(DE.d,1)'; Par_info.max(1:end) = DE.exponents(end); 
Par_info.boundhandling = 'reflect';        % Explicit boundary handling
% Create initial population
pop = nan (DE.N, DE.d);
for ii = 1:DE.N
    % create a random population with (true) replacement
     pop(ii,:) = randsample(DE.exponents,DE.d,'true');
end
% DE Algorithm
[Output_MODE] = diff_evolution_MO (Y_cal, X_cal, DE, pop, Par_info);

% This function executes a differential evolution nonlinear regression
% Author: Guilherme Jose Cunha Gomes
% date: 08 / MAY / 2015
function [Output_DE] = diff_evolution_MO (Y_cal, X_cal, DE, pop, Par_info)
%% Pre-processing
 
% Prelocate memory for objective function matrix
store_OF = nan(DE.last,2);

% Prelocate memory for storage of best parameters
   
if strcmp(DE.bias, 'Yes')% Bias considered 
    store_PAR = nan(DE.last, DE.m + 1);
else %No bias
    store_PAR = nan(DE.last, DE.m);
end

% Prelocate memory for model structure of the best individual
store_POP = nan(DE.last, DE.d);

% Prelocate memory for storage of output of the best individual 
store_y = nan(DE.last, DE.n); 

% Prelocate memory for storage of output of the correlations between
% dependent and indepent variables
store_r = nan(DE.last, size(X_cal,2)); 

% Prelocate memory to store statistics for R2, RMSE, MAE, r, E_rel
store_stat = nan(DE.last, 5); 

% setup compromise programing
C = strings(1,size(store_OF,2))';  % create matrix of strings
C(1:end) = 'min';                  % minimize all terms
a = ones(1,numel(C))';             % equal weights


%% Start EPR

% start counting
%tic

% waitbar
%hw = waitbar(0,'Running MODE...');

    % Start evolution
    for ib = 1:DE.last
        % update wait bar
         %waitbar(ib/DE.last);
        % Evaluate population
        [y, theta] = evolutionary_model(pop, DE, X_cal, Y_cal);
        % Differential evolution tool 
        [pop, OF, y, theta] = DE_algorithm (pop, DE, Par_info, X_cal, Y_cal, y, theta, store_POP, ib);
        % compromise programing to store best solutions of objective functions
        [Out_MO] = compromise (OF', C, a);
        % Store "best" simulation
        store_y(ib,:) = y(:,Out_MO.Y);
        % Store "best" values of objective function
        store_OF(ib,:) = OF(Out_MO.Y,:);
        % Store "best" values of EPR parameters
        store_PAR(ib,:) = theta(:,Out_MO.Y);
        % Store model structure of the best individual 
        store_POP(ib,:) = pop(Out_MO.Y,:);
        % store correlation between Y-output and X-variables
        for jj = 1:size(X_cal,2)
           [Out_stat] = statistic_model (X_cal(:,jj), y(:,Out_MO.Y));
            store_r(ib,jj) = Out_stat.r;
        end
        % compute and store statistics of the best simulation
        [Out_stat] = statistic_model (Y_cal, y(:,Out_MO.Y));
        store_stat(ib,1) = Out_stat.rsq;    % R2
        store_stat(ib,2) = Out_stat.RMSE;   % RMSE
        store_stat(ib,3) = Out_stat.MAE;    % MAE
        store_stat(ib,4) = Out_stat.E_rel;  % E_rel
        store_stat(ib,5) = Out_stat.r;      % r

    end
    % close waitbar
   % close(hw);
   % toc

    %% Post processing

    % Select best population
    bestpop = store_POP(end,:);

    % Select best output
    besty = store_y(end,:)';

    % Select best parameters
    best_par = store_PAR(end,:);

    % compute statistics
    [Out_stat] = statistic_model (Y_cal, besty); 

    % Create a structured array for the final solution
    Output_DE.bestpop = (reshape(bestpop, DE.k, DE.m))'; % best population
    Output_DE.besty = besty;                             % best simulation
    Output_DE.best_par = best_par;                       % best parameters  
    Output_DE.RMSE = Out_stat.RMSE;                      % RMSE  
    Output_DE.MAE = Out_stat.MAE;                        % MAE
    Output_DE.R2 = Out_stat.R2;                          % R2
    Output_DE.r = Out_stat.r;                            % r
    Output_DE.sse_min = Out_stat.SSE;                    % minimum SSE
    Output_DE.E_rel = Out_stat.E_rel;                    % E_rel
    Output_DE.OF = store_OF;                             % Objective functions
    Output_DE.SEE = store_OF(:,1);                       % store SSE
    Output_DE.stats = store_stat;                        % store statistics
                     
end

function [pop, OF, y, theta] = DE_algorithm (pop, DE, Par_info, X, Y, y, theta, store_POP, ib)

% Differential evolution tool %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % algorithmic variables
     N = DE.N; delta = DE.delta; R = DE.R; gamma = DE.gamma; cr = DE.cr;

     % Randomly permute [1,...,N?1] N times
     [~,draw] = sort(rand(N-1,N));
          
     % create offspring from parent population
     z = zeros(size(pop));
     for i = 1:N
         
         r1 = R(i,draw(1:delta,i));         % Derive vector r1
         r2 = R(i,draw(delta+1:2*delta,i)); % Derive vector r2
         % generate new population
         z (i,:) = pop(i,:) + gamma * (pop(r1,:) - pop(r2,:)); 
        
         % now deal with boundaries ('bound')
          [z] = Boundary_handling(z,Par_info);
         
         % Modify only dimensions of "pop" if u > crossover probability, CR
         A = find(rand(1,DE.d) < cr);
         z (i,A) = pop(i,A);
         
         % check if exist any EPR term without a value (optional). With
         % this function, one can make sure that at least on term of the
         % EPR will recieve a variable and at least one variable will have
         % a non-zero exponent. The presence of at least one zero in EX 
         % ensures the ability to exclude some of the inputs and/
         % or input combinations from the regression equation. If you are
         % sure that the variables used must be ib the model, you can use 
         % the function below. Otherwise, comment this function. 
         [z] = check_variable(ib,z,i,DE,store_POP);

     end
    
    % Evaluate model for offspring
    [y_off, theta_off] = evolutionary_model(z, DE, X, Y);
    
    % create vectors to store objective functions of parents
    SSE = zeros(1,size(y,2));  % SSE
    for kk = 1:size(y,2)
        [Out_stat] = statistic_model (y(:,kk), Y);
        SSE(kk) = Out_stat.SSE;
    end
    
    % minimization of the number of inputs
    comb = ones(size(pop));    % create matrix to count number of inputs
    id = find(pop==0);         % find no input values
    comb(id) = 0;              % update comb matrix
    r = sum(comb,2);           % now sum the combinations as objective target

    % Concatenate objective function of parents 
    Xp = [SSE' r];
    
    % create vectors to store objective functions of offspring
    SSE = zeros(1,size(y,2));  % SSE
    for kk = 1:size(y,2)
        [Out_stat] = statistic_model (y_off(:,kk), Y);
        SSE(kk) = Out_stat.SSE;
    end
    % minimization of the number of inputs
    comb = ones(size(z));    % create matrix to count number of inputs
    id = find(z==0);         % find no input values
    comb(id) = 0;              % update comb matrix
    r = sum(comb,2);           % now sum the combinations as objective target 
    
    % Concatenate objective function of parents 
    Xo = [SSE' r];
    
    % Now we put all objective functions in one matrix to be evaluated
      X = [Xp; Xo]; 
 
    % now rank objective functions of parents and offspring
    [ rank , CrowdDist ] = pareto_ranking ( X );
    
    % reshape rank in order to compare parents and offspring ranks
    new_rank = reshape(rank,DE.N,2);
    
    % evaluate rank (find ranks in offspring that are lower than parents)
    dx = find(new_rank(:,2) < new_rank(:,1));
    
    % now change new population
    pop(dx,:) = z(dx,:);
    
    % correct eventual population near zero value
    idx = find(abs(pop)<0.01);
    pop(idx) = 0;
    
    % update objective function
    Xp(dx,:) = Xo(dx,:); OF = Xp;
    
    % update EPR parameters
    theta(:,dx) = theta_off(:,dx);
    
    % update model response
    y(:,dx) = y_off(:,dx);
    
end

    function [z] = check_variable(ib,z,i,DE,store_POP)
    % make sure every EPR term will have at least one variable with a
    % non zero exponent
        if ib>1 && any(sum(abs(reshape(z(i,:),DE.k,DE.m)'),2) == 0)
            z (i,:) = store_POP(ib-1,:);    % convert to the best simulation
            b = randi([1 numel(z(i,:))]);   % perturb the best simulation
            expon = DE.exponents;           % check exponents to 
            [minValue, closestIndex] = min(abs(z(i,b) - expon));
            expon(closestIndex) = [];
            z(i,b) = randsample(expon,1);
        end
        % make sure that each variable will have one exponent (at least)
        if ib>1 && numel(z(i,:))>DE.k && any(sum(abs(reshape(z(i,:),DE.k,DE.m)'),1) == 0)
            % store best simulation from previous generation
            ne = store_POP(ib-1,:); 
            % identify a pool of vectors to sample and remove zero
            % values
            exp_vec = DE.exponents; exp_vec((exp_vec==0))=[];
            % make sure that each variable will have one exponent (at least)
            if any(sum(abs(reshape(ne,DE.k,DE.m)'),1) == 0)
                % find columns with all zeros'
                cols_zeros = find(all(ne==0));
                if numel(cols_zeros)==1
                % select a row to change the exponent
                sel_row = randi([1 numel(ne(:,cols_zeros))]);
                % update population
                ne(sel_row,cols_zeros)=randsample(exp_vec,1,'true');
                z (i,:) = ne;
                else
                    for jj = 1:numel(cols_zeros)
                        % select a row to change the exponent
                        sel_row = randi([1 numel(ne(:,jj))]);
                        % update population
                        ne(sel_row,jj)=randsample(exp_vec,1,'true');
                        z (i,:) = ne;
                    end
                end
            else % change olny one exponent of the best simulation
                % select an exponent to change 
                sel_exp = randsample(1:numel(ne),1,'true');
                % update population
                ne(sel_exp)=randsample(exp_vec,1,'true');
                z (i,:) = ne;
            end
                
        end    
    end

end    