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
     
     % check if exist any EPR term without a value
     for i = 1:N
         if any(sum(abs(reshape(newpop(i,:),GA.k,GA.m)'),2) == 0)
             newpop (i,:) = pop(i,:);
         end
         if any(sum(abs(reshape(newpop(i,:),GA.k,GA.m)'),1) == 0)
            newpop (i,:) = pop(i,:);
        end
     end
     % Update new population;
     pop = newpop;
     
end