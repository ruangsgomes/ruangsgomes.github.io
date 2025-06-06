function [ rank , CrowdDist ] = pareto_ranking ( X );
% Calculate the Pareto rank of each solution of X
%
% Input: matrix of size N by m, where:
%               N is the number of individuals, 
%               m are their objective function values (m > 1)
%
% Output:
% 
% "rank"        --> Pareto rank
% "CrowdDist"   --> Crowding distance; as large as possible
%
% Written by Jasper A. Vrugt
%
% Example:
% 
% X = rand(100,2); % --> Create 100 individuals with 2 objective function values
% [ rank , CrowdDist ] = pareto_rank ( X );
%
% If code is not working --> then recompile the paretofront.c code using
% mex paretofront.c
% This only needs to be done once. 

% Determine the size of X
[Nx Ny] = size ( X ); 

% Initialize the vector rank with zeros
rank = NaN ( Nx , 1 ); 

% Initialize Crowding distance
CrowdDist = NaN ( Nx , 1 );

% Initialize rank of the front
rank_of_front = 1;

% Check whether all solutions have been ranked?
idx = find ( isnan ( rank ) );

% Now a while statement --> as long as not all solutions have been ranked
while isempty(idx) == 0,
    
    % Now evaluate these points 
    front = paretofront ( X ( idx , 1 : Ny ) );

    % Which values are true
    ii = idx ( front == 1 ); N = numel(ii); 
    
    % Update rank
    rank ( ii ) = rank_of_front;
    
    % Which solutions of X are not ranked yet?
    idx = find ( isnan ( rank ) );

    % Now select objective functions 
    F = X ( ii , 1 : Ny );
    
    % Initialize distance
    D = zeros ( N , 1);
    
    % Now loop over objectives
    for j = 1 : Ny,
        
        % Sort the jth function values
        [ Fs , s_idx ] = sort ( F ( : , j ) );
    
        % Assing extreme values to end
        D ( s_idx([1 N]) ) = D( s_idx([1 N]) ) + Inf;
        
        % Now compute difference of two points
        D ( s_idx( 2 : N - 1 ) ) =  D ( s_idx( 2 : N - 1) ) + ( Fs ( 3 : N ) - Fs ( 1 : N - 2 ) );

    end;
    
    % Store crowding distance
    CrowdDist ( ii , 1 ) = D;
    
    % Update rank of the front
    rank_of_front = rank_of_front + 1;
    
end;