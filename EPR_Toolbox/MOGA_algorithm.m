function [Output_MOGA] = MOGA_algorithm (GA, Y_cal, X_cal)

%% MOGA - Multi-objective genetic algorithm

% dimension of the problem
GA.d = GA.m*GA.k;

% Create initial population
pop = nan (GA.N, GA.d);
for ix = 1:GA.N
    % create a random population with (true) replacement
    pop(ix,:) = randsample(GA.exponents,GA.d,'true');
end

% This function executes a Multi-objective genetic algorithm EPR
[Output_MOGA] = gen_alg_MO (Y_cal, X_cal, GA, pop);


function [Output_MOGA] = gen_alg_MO (Y_cal, X_cal, GA, pop)
% Author: Guilherme Jose Cunha Gomes
% date: 08 / MAY / 2015
%% initial computation


% Prelocate memory for objective function matrix
store_OF = nan(GA.last,2);
% SSE plus correlation coefficients between independent variables and
% dependent Y variable
%store_OF = nan(DE.last,size(X_cal,2)+1);

% Prelocate memory for storage of best parameters
if strcmp(GA.bias, 'Yes')% Bias considered 
    store_PAR = nan(GA.last, GA.m + 1);
else %No bias
    store_PAR = nan(GA.last, GA.m);
end  
%store_PAR = nan(GA.last, GA.m + 1);

% Prelocate memory for model structure of the best individual
store_POP = nan(GA.last, GA.d);

% Prelocate memory for storage of output of the best individual 
store_y = nan(GA.last, GA.n); 

% Prelocate memory for storage of output of the correlations between
% dependent and indepent variables
store_r = nan(GA.last, size(X_cal,2)); 

% Prelocate memory to store statistics for R2, RMSE, MAE, r, E_rel
store_stat = nan(GA.last, 5); 
% setup compromise programing
C = strings(1,size(store_OF,2))';  % create matrix of strings
C(1:end) = 'min';                  % minimize all terms
a = ones(1,numel(C))';             % equal weights


%% genetic algorithm evolution

% start counting
%tic
% waitbar
%hw = waitbar(0,'Running MOGA...');
% Start evolution
for ib = 1:GA.last
    % update wait bar
    % waitbar(ib/GA.last);
    % Evaluate population
    [y, theta] = evolutionary_model(pop, GA, X_cal, Y_cal);
    % genetic algorithm tool 
    [pop, OF, y, theta] = gen_algorithm_MO (pop, y, GA, X_cal, Y_cal, theta, store_POP, ib);
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

%close(hw);
%toc

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
Output_MOGA.bestpop = (reshape(bestpop, GA.k, GA.m))'; % best population
Output_MOGA.besty = besty;                             % best simulation
Output_MOGA.best_par = best_par;                       % best parameters
Output_MOGA.RMSE = Out_stat.RMSE;                      % RMSE 
Output_MOGA.MAE = Out_stat.MAE;                        % MAE 
Output_MOGA.R2 = Out_stat.R2;                          % R2
Output_MOGA.r = Out_stat.r;                            % r
Output_MOGA.sse_min = Out_stat.SSE;                    % minimum SSE
Output_MOGA.E_rel = Out_stat.E_rel;                    % E_rel
Output_MOGA.OF = store_OF;                             % Objective functions
Output_MOGA.SEE = store_OF(:,1);                       % store SSE
Output_MOGA.stats = store_stat;                        % store statistics

end


% % This function is a global optimization tool using Genetic Algorithm (GA)
function [pop, OF, y, theta] = gen_algorithm_MO (pop, y, GA, X_cal, Y_cal, theta, store_POP, ib)
% Author: Guilherme Gomes 

% GA settings:
%     selection: Normalized geometric ranking selection
%     crossover: Multiple point crossover
%     mutation: Single point mutation

% setup algorithmic parameters 
m = GA.m; k = GA.k; cpr = GA.cpr; mutrate = GA.mutrate; N = GA.N;

% Evaluate parents %%%%%%%%%%%%%%%%%%%%%%%%%
% create vectors to store objective functions of parents
SSE = zeros(1,size(y,2));  % SSE
for kk = 1:size(y,2)
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

% Now sort SSE
[SSE_sort,ind]=sort(SSE);

% Selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
% Normalized geometric ranking selection
% P[selecting the ith individual] = q^t(1-q)^(r-1)
     
% Compute probability
prob = (1./SSE)/sum(1./SSE);
% find maximum individual probability
q = max(prob);
% normalize
ql = q/(1-((1-q)^N));
% create auxiliar variables
P = zeros(N,1); PP = zeros(N,1);
% now normalize all population values. The first will have the higher
% probability, the last will have the worst.
for i = 1:N
    P(i,1) = ql*(1-q)^(i-1);
end
% we attribute to PP the index of fitness of the objective function
PP(ind(1:N))=P(1:N);
% After this point PP is vector of probabilities, and we take the 
% cumulative sum
roulette = cumsum(PP);
% choose N randon numbers of the uniform distribution
rn = rand(N,1);
% selection by roulette wheel
idx = zeros(N,1);
for ii=1:N;
    idx(ii,:) = find(roulette > rn(ii), 1, 'first');
end

% update new population
newpop = pop(idx,:);
     
% Crossover %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first find points to perform crossover
rcross = rand(N,1);
selcross = find(rcross<cpr);
% changing chromosomes selected for even number
if mod(numel(selcross),2)==1;
    unselect = randi(numel(selcross));
    selcross(unselect,:) = [];
end
% multiple point crossover
for i = 1:2:numel(selcross)
    % select two chromosomes randomly to perform crossover
    yy = randsample(numel(selcross),2);
    % identify father and mother
    father = newpop(selcross(yy(1)),:);
    mother = newpop(selcross(yy(2)),:);
    % select number of crossover points
    mp = randi((k*m)-1);
    % choosing the points of crossover
    cpt = randsample(k*m-1,mp); cpt=sort(cpt); % crossover points
    for j = 1:mp-1
        newpop(selcross(yy(1)),cpt(j):cpt(j+1))=mother(cpt(j):cpt(j+1));
        newpop(selcross(yy(2)),cpt(j):cpt(j+1))=father(cpt(j):cpt(j+1));
    end
    selcross(yy,:) = [];
end
     
% Mutation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select points to perform mutation 
rmut = rand(N,m*k);
[rows,cols] = find(rmut<mutrate);
for i = 1:numel(rows)
    newpop(rows(i),cols(i))=randsample(GA.exponents,1);
end

% check if exist any EPR term without a value (optional). With
% this function, one can make sure that at least on term of the
% EPR will recieve a variable and at least one variable will have
% a non-zero exponent. The presence of at least one zero in EX 
% ensures the ability to exclude some of the inputs and/
% or input combinations from the regression equation. If you are
% sure that the variables used must be ib the model, you can use 
% the function below. Otherwise, comment this function. 
for i = 1:N
    [newpop] = check_variable(ib,newpop,i,GA,store_POP);
end

% end gen algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Fitness function of Ofspring %%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate model for offspring
[y_off, theta_off] = evolutionary_model(newpop, GA, X_cal, Y_cal);
    
% create vectors to store objective functions of parents
SSE = zeros(1,size(y,2));  % SSE
for kk = 1:size(y,2)
    [Out_stat] = statistic_model (y_off(:,kk), Y_cal);
    SSE(kk) = Out_stat.SSE;
end
    
% minimization of the number of inputs
comb = ones(size(newpop));    % create matrix to count number of inputs
id = find(newpop==0);         % find no input values
comb(id) = 0;              % update comb matrix
r = sum(comb,2);           % now sum the combinations as objective target

% Concatenate objective function of parents 
Xo = [SSE' r];  
 
% Now we put all objective functions in one matrix to be evaluated
X = [Xp; Xo];
    
% now rank objective functions of parents and offspring
[ rank , CrowdDist ] = pareto_ranking ( X );
    
% reshape rank in order to compare parents and offspring ranks
new_rank = reshape(rank,GA.N,2);
    
% evaluate rank (find ranks in offspring that are lower than parents)
dx = find(new_rank(:,2) < new_rank(:,1));
    
% now change new population
pop(dx,:) = newpop(dx,:);
% correct evaentual population near zero value
idx = find(abs(pop)<0.01);
pop(idx) = 0;
    
% update objective function
Xp(dx,:) = Xo(dx,:); OF = Xp; 
% update EPR parameters
theta(:,dx) = theta_off(:,dx);
% update model response
y(:,dx) = y_off(:,dx);     
end    


function [newpop] = check_variable(ib,newpop,i,GA,store_POP)
    % check if exist any EPR term without a value
    if ib>1 && any(sum(abs(reshape(newpop(i,:),GA.k,GA.m)'),2) == 0)
        newpop (i,:) = store_POP(ib-1,:);    % convert to the best simulation
            b = randi([1 numel(newpop(i,:))]); % perturb the best simulation
            expon = GA.exponents;         % check exponents to 
            [minValue, closestIndex] = min(abs(newpop(i,b) - expon));
            expon(closestIndex) = [];
            newpop(i,b) = randsample(expon,1);
          %  best simulation store_POP(ib,:);
    end
    % make sure that each variable will have one exponent (at least)
    if ib>1 && numel(newpop(i,:))>GA.k && any(sum(abs(reshape(newpop(i,:),GA.k,GA.m)'),1) == 0)
            % store best simulation from previous generation
            ne = store_POP(ib-1,:); 
            % identify a pool of vectors to sample and remove zero
            % values
            exp_vec = GA.exponents; exp_vec((exp_vec==0))=[];
            % make sure that each variable will have one exponent (at least)
            if any(sum(abs(reshape(ne,GA.k,GA.m)'),1) == 0)
                % find columns with all zeros'
                cols_zeros = find(all(ne==0));
                if numel(cols_zeros)==1
                % select a row to change the exponent
                sel_row = randi([1 numel(ne(:,cols_zeros))]);
                % update population
                ne(sel_row,cols_zeros)=randsample(exp_vec,1,'true');
                newpop (i,:) = ne;
                else
                    for jj = 1:numel(cols_zeros)
                        % select a row to change the exponent
                        sel_row = randi([1 numel(ne(:,jj))]);
                        % update population
                        ne(sel_row,jj)=randsample(exp_vec,1,'true');
                        newpop (i,:) = ne;
                    end
                end
            else % change olny one exponent of the best simulation
                % select an exponent to change 
                sel_exp = randsample(1:numel(ne),1,'true');
                % update population
                ne(sel_exp)=randsample(exp_vec,1,'true');
                newpop (i,:) = ne;
            end
     end
end
     
end % end MOGA