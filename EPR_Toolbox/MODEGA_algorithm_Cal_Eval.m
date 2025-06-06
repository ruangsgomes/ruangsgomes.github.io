function [Output_MODEGA] = MODEGA_algorithm_Cal_Eval (DE, GA, Y_cal, X_cal,Y_eval, X_eval)
%function [Output_MODEGA] = MODEGA_algorithm (DE, GA, Y_cal, X_cal)
%% MODEGA - Multi-objective Differential evolution algorithm

% dimension of the problem
DE.d = DE.m*DE.k;  
% Parameter information
Par_info.min = zeros(DE.d,1)'; Par_info.min(1:end) = DE.exponents(1);  
Par_info.max = zeros(DE.d,1)'; Par_info.max(1:end) = DE.exponents(end); 
Par_info.boundhandling = 'reflect';        % Explicit boundary handling

% MODEGA Algorithm
[Output_MODEGA] = DE_GA_regression (Y_cal, X_cal,Y_eval, X_eval, DE, GA, Par_info);
% This function executes a dual search based EPR modeling approach
% Author: Guilherme Jose Cunha Gomes
% date: 08 / MAY / 2015
end

function [Output_MODEGA] = DE_GA_regression (Y_cal, X_cal,Y_eval, X_eval, DE, GA, Par_info)
%function [Output_MODEGA] = DE_GA_regression (Y_cal, X_cal, DE, GA, Par_info)
%% initial settings
% create a vector to store population size for both GA and DE
% during different generations
N_alg = nan(DE.last,2);

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

% Prelocate memory for storage of output of the best individual 
store_y_eval = nan(DE.last, size(X_eval,1)); 

% Prelocate memory for storage of output of the correlations between
% dependent and indepent variables
store_r = nan(DE.last, size(X_cal,2)); 

% Prelocate memory to store statistics for R2, RMSE, MAE, r, E_rel
store_stat = nan(DE.last, 5); 

% Prelocate memory to store evaluation statistics for R2, RMSE, MAE, r, E_rel
store_stat_eval = nan(DE.last, 5); 

% setup compromise programing
C = strings(1,size(store_OF,2))';  % create matrix of strings
C(1:end) = 'min';                  % minimize all terms
a = ones(1,numel(C))';             % equal weights


%% initial population 
% R-matrix for index     
for i = 1:DE.N, DE.R(i,1:DE.N-1) = setdiff(1:DE.N,i); end 
% Create initial population
pop = nan (DE.N, DE.d);
for i = 1:DE.N
    % create a random population with (true) replacement
    pop(i,:) = randsample(DE.exponents,DE.d,'true');
end

%% Run Evolutionary model
[y, theta] = evolutionary_model(pop, DE, X_cal, Y_cal);

%% create vector to store objective function
SSE = zeros(1,DE.N);  % SSE
% Evaluate fitness
for kk = 1:DE.N
    [Out_stat] = statistic_model (y(:,kk), Y_cal);
    SSE(kk) = Out_stat.SSE;
end
% minimization of the number of inputs
comb = ones(size(pop));    % create matrix to count number of inputs
id = find(pop==0);         % find no input values
comb(id) = 0;              % update comb matrix
r = sum(comb,2);           % now sum the combinations as objective target

% Concatenate objective function of parents 
Xp = [SSE' r];

%% Create first offspring generation
% Differential evolution tool %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = round(DE.N/2); % 50% percent 
delta = DE.delta; R = DE.R; gamma = DE.gamma; cr = DE.cr;
% Randomly permute [1,...,N?1] N times
[~,draw] = sort(rand(N-1,N));
% create offspring from parent population
z = zeros(N,numel(pop(1,:)));
for i = 1:N
    r1 = R(i,draw(1:delta,i));         % Derive vector r1
    r2 = R(i,draw(delta+1:2*delta,i)); % Derive vector r2
    z (i,:) = pop(i,:) + gamma * (pop(r1,:) - pop(r2,:)); 
    % now deal with boundaries ('bound')
    [z] = Boundary_handling(z,Par_info);
    % Modify only dimensions of "pop" if u > crossover probability, CR
    A = find(rand(1,DE.d) < cr);
    z (i,A) = pop(i,A);
    % check if exist any EPR term without a value
    if any(sum(abs(reshape(z(i,:),DE.k,DE.m)'),2) == 0)
       z (i,:) = pop(i,:); 
    end
    if any(sum(abs(reshape(z(i,:),DE.k,DE.m)'),1) == 0)
       z (i,:) = pop(i,:); 
    end
end

% Now sort SSE
pop_ga = pop(N+1:end,:);   % select population for offspring genetation
SSE_ga = SSE(N+1:end);     % select objective function
[~,ind]=sort(SSE_ga);
GA.N = size(pop_ga,1);
% execute genetic operator with 50% of the population
[pop_ga_off] = gen_algorithm (pop_ga, SSE_ga, ind, GA);
%[pop_ga_off] = gen_algorithm_MO (pop_ga, SSE_ga, ind, GA, store_POP, ib);
% concanetate population of the first generation
pop_off = [z; pop_ga_off];

% Run Evolutionary model
[y_off, theta_off] = evolutionary_model(pop_off, DE, X_cal, Y_cal);

% create vector to store objective function
SSE = zeros(1,DE.N);  % SSE
% Evaluate fitness
for kk = 1:DE.N
    [Out_stat] = statistic_model (y_off(:,kk), Y_cal);
    SSE(kk) = Out_stat.SSE;
end
  
% minimization of the number of inputs
comb = ones(size(pop_off));    % create matrix to count number of inputs
id = find(pop_off==0);         % find no input values
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
pop(dx,:) = pop_off(dx,:);
% correct evaentual population near zero value
idx = find(abs(pop)<0.01);
pop(idx) = 0;
    
% update objective function
Xp(dx,:) = Xo(dx,:); OF = Xp; 
% update EPR parameters
theta(:,dx) = theta_off(:,dx);
% update model response
y(:,dx) = y_off(:,dx);

%-------------- Self adaptive offspring creation
Nde = N; Nga = N;
N_alg(1,:) = [Nde Nga];  % store population size
N_de_new = numel(find(dx<=N));

if (numel(find(dx<=N))/N) < 0.2
    Nde = round(DE.N * 0.2);  % population size for DE
    Nga = DE.N - Nde;  % population size for GA
else
    Nde = round(DE.N*(N_de_new)/N_alg(1,1));  % population size for DE
    Nga = DE.N - Nde;                   % population size for GA
end


% compromise programing to store best solutions of objective functions
[Out_MO] = compromise (OF', C, a);
% first simulation index
ib = 1;
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
store_stat(ib,6) = Out_stat.PBIAS;  % PBIAS

% Use best EPR structure to compute evaluation statistics
Output_eval.bestpop = store_POP(ib,:);
Output_eval.bestpop = (reshape(Output_eval.bestpop, DE.k, DE.m))';

Output_eval.best_par = store_PAR(ib,:);

[y_eval] = evaluation_DE_regression (DE, Output_eval, X_eval);

% compute and store statistics for evaluation dataset
[Out_stat_eval] = statistic_model (Y_eval, y_eval);
store_stat_eval(ib,1) = Out_stat_eval.rsq;    % R2
store_stat_eval(ib,2) = Out_stat_eval.RMSE;   % RMSE
store_stat_eval(ib,3) = Out_stat_eval.MAE;    % MAE
store_stat_eval(ib,4) = Out_stat_eval.E_rel;  % E_rel
store_stat_eval(ib,5) = Out_stat_eval.r;      % r
store_stat_eval(ib,6) = Out_stat_eval.PBIAS;  % PBIAS
    
%% Algorithm evolution
%hw = waitbar(0,'Running MODEGA EPR...');
% Start evolution
for ib = 2:DE.last
    % update wait bar
    % waitbar(ib/DE.last);
    % Create first offspring generation
    % Differential evolution tool %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = N_alg(ib-1,1); % size population for DE 
    % Randomly permute [1,...,N?1] N times
    [~,draw] = sort(rand(N-1,N));
    % create offspring from parent population
    z = zeros(N,numel(pop(1,:)));
    for i = 1:N
        r1 = R(i,draw(1:delta,i));         % Derive vector r1
        r2 = R(i,draw(delta+1:2*delta,i)); % Derive vector r2
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

    % Now sort SSE
    pop_ga = pop(N+1:end,:);   % select population for offspring genetation
    SSE_ga = OF(N+1:end,1);    % select SSE
    [~,ind]=sort(SSE_ga);
    GA.N = size(pop_ga,1);
    %[pop_ga_off] = gen_algorithm_MO (pop_ga, SSE_ga, ind, GA);
    [pop_ga_off] = gen_algorithm (pop_ga, SSE_ga, ind, GA);
    %[pop_ga_off] = gen_algorithm_MO (pop_ga, SSE_ga, ind, GA, store_POP, ib);
    % concanetate population of the first generation
    pop_off = [z; pop_ga_off];

    % Run Evolutionary model
    [y_off, theta_off] = evolutionary_model(pop_off, DE, X_cal, Y_cal);

    % create vector to store objective function
    SSE = zeros(1,DE.N);  % SSE
    % Evaluate fitness
    for kk = 1:DE.N
        [Out_stat] = statistic_model (y_off(:,kk), Y_cal);
        SSE(kk) = Out_stat.SSE;
    end
   
    % minimization of the number of inputs
    comb = ones(size(pop_off));    % create matrix to count number of inputs
    id = find(pop_off==0);         % find no input values
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
    pop(dx,:) = pop_off(dx,:);
    % correct evaentual population near zero value
    idx = find(abs(pop)<0.01);
    pop(idx) = 0;

    % update objective function
    Xp(dx,:) = Xo(dx,:); OF = Xp; 
    % update EPR parameters
    theta(:,dx) = theta_off(:,dx);
    % update model response
    y(:,dx) = y_off(:,dx);

    % Self adaptive offspring creation
    N_de_new = numel(find(dx<=N));
    if sum(dx)==0
        Nde = N; Nga = DE.N - Nde;
    elseif (numel(find(dx<=N))/numel(dx)) < 0.2 
        Nde = round(DE.N * 0.2);  % population size for DE
        Nga = DE.N - Nde;  % population size for GA
    elseif (numel(find(dx<=N))/numel(dx)) > 0.8
        Nde = round(DE.N * 0.8);  % population size for DE
        Nga = DE.N - Nde;  % population size for GA
    else
        Nde = round(DE.N*(N_de_new)/numel(dx));  % population size for DE
        Nga = DE.N - Nde;                   % population size for GA
    end
    N_alg(ib,:) = [Nde Nga];  % store population size

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
    store_stat(ib,6) = Out_stat.PBIAS;  % PBIAS
    
    % Use best EPR structure to compute evaluation statistics
    Output_eval.bestpop = store_POP(ib,:);
    Output_eval.bestpop = (reshape(Output_eval.bestpop, DE.k, DE.m))';
    Output_eval.best_par = store_PAR(ib,:);

    [y_eval] = evaluation_DE_regression (DE, Output_eval, X_eval);

    % compute and store statistics for evaluation dataset
    [Out_stat_eval] = statistic_model (Y_eval, y_eval);
    store_stat_eval(ib,1) = Out_stat_eval.rsq;    % R2
    store_stat_eval(ib,2) = Out_stat_eval.RMSE;   % RMSE
    store_stat_eval(ib,3) = Out_stat_eval.MAE;    % MAE
    store_stat_eval(ib,4) = Out_stat_eval.E_rel;  % E_rel
    store_stat_eval(ib,5) = Out_stat_eval.r;      % r
    store_stat_eval(ib,6) = Out_stat_eval.PBIAS;  % PBIAS
end

%close(hw);
%% Post processing

% Select best population
bestpop = store_POP(end,:);

% Select best output
besty = store_y(end,:)';

% Select best parameters
best_par = store_PAR(end,:);

% compute statistics
[Out_stat] = statistic_model (Y_cal, besty); 

% Create a solution structured array
    % Create a structured array for the final solution
    Output_MODEGA.bestpop = (reshape(bestpop, DE.k, DE.m))'; % best population
    Output_MODEGA.besty = besty;                             % best simulation
    Output_MODEGA.best_par = best_par;                       % best parameters  
    Output_MODEGA.RMSE = Out_stat.RMSE;                      % RMSE  
    Output_MODEGA.MAE = Out_stat.MAE;                        % MAE
    Output_MODEGA.R2 = Out_stat.R2;                          % R2
    Output_MODEGA.r = Out_stat.r;                            % r
    Output_MODEGA.sse_min = Out_stat.SSE;                    % minimum SSE
    Output_MODEGA.E_rel = Out_stat.E_rel;                    % E_rel
    Output_MODEGA.PBIAS = Out_stat.PBIAS;                    % PBIAS
    Output_MODEGA.OF = store_OF;                             % Objective functions
    Output_MODEGA.SEE = store_OF(:,1);                       % store SSE
    Output_MODEGA.stats = store_stat;                        % store statistics
    Output_MODEGA.stats_eval = store_stat_eval;                        % store statistics
    Output_MODEGA.N_alg = N_alg;                             % number of offsprings
end % MODEGA function



function [z] = check_variable(ib,z,i,DE,store_POP)
    % make sure every EPR term will have at least one variable with a
    % non zero exponent (optional)
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