% This function executes a genetic algorithm nonlinear regression

% Author: Guilherme Jose Cunha Gomes

% date: 07 / AUG / 2019

function [Output_GA] = gen_alg_regression (pop,GA,X_cal,Y_cal)

    % loop throught all generations
    for ib = 1:GA.last

        % Evaluate model
        [y, theta] = evolutionary_model(pop, GA, X_cal, Y_cal);

        % Now create a matrix of residuals
        Res = y - Y_cal; 

        % compute objective function (sum of squared errors : SSE)
        SSE = sum(Res.^2)./GA.n;

        % Also compute Coeffient of Determination
        CoD = 1 - ((SSE./((sum(y,1) - mean(y,1)).^2))*GA.n);

        % Now sort SSE
        [SSE_sort,ind]=sort(SSE);

        % Store objective function of the best individual 
        store_OF(ib,:) = [SSE_sort(1) CoD(ind(1))];

        % Store parameters of the best individual 
        store_PAR(ib,:) = theta(:,ind(1));

        % Store model structure of the best individual 
        store_POP(ib,:) = pop(ind(1),:);

        % Store output of the best individual 
        store_y(ib,:) = y(:,ind(1));


        % Check convergence --> single objective %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %     % Stop running if SSE is zero 
    %      if SSE_sort(1)<=1e-1
    %          break
    %      end


        % Genetic algorithm tool %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [pop] = gen_algorithm (pop, SSE, ind, GA);
    end

toc

% Post processing

% Find min SSE and max CoD
[ro, co] = min(store_OF(:,1));

% Select best population
bestpop = store_POP(co,:);

% Select best output
besty = store_y(co,:)';

% Select best parameters
best_par = store_PAR(co,:);

% compute statistics
[Out_stat] = statistic_model (Y_cal, besty); 

% Create a solution structured array
%Output_GA.bestpop = reshape(bestpop, GA.m, GA.k);
Output_GA.bestpop =(reshape(bestpop, GA.k, GA.m))';
Output_GA.besty = besty;
Output_GA.best_par = best_par;
Output_GA.rsq = Out_stat.rsq;
Output_GA.RMSE = Out_stat.RMSE;
Output_GA.R2 = Out_stat.R2;
Output_GA.r = Out_stat.r;
Output_GA.MAE = Out_stat.MAE;
Output_GA.sse_min = Out_stat.SSE;
Output_GA.E_rel = Out_stat.E_rel;
Output_GA.CoD = store_OF(:,2);
Output_GA.SEE = store_OF(:,1);


    % store correlation between Y-output and X-variables
    for jj = 1:size(X_cal,2)
        [Out_stat] = statistic_model (X_cal(:,jj), besty);
        store_r(1,jj) = Out_stat.r;
    end
Output_GA.correl = store_r;    
end

function [pop] = gen_algorithm (pop, SSE, ind, GA)

m = GA.m; k = GA.k; cpr = GA.cpr; mutrate = GA.mutrate; N = GA.N;
% This function is a global optimization using Genetic Algorithm (GA)

% Author: Guilherme Gomes 

% GA settings:
%     selection: Normalized geometric ranking selection
%     crossover: Multiple point crossover
%     mutation: Single point mutation


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

     % New population
     newpop = pop(idx,:);
     
% Crossover %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
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
     % 
     rmut = rand(N,m*k);
     [rows,cols] = find(rmut<mutrate);
     for i = 1:numel(rows)
         newpop(rows(i),cols(i))=randi([GA.exponents(1) GA.exponents(end)]);
     end
      
     % Update new population;
     pop = newpop;
     
end