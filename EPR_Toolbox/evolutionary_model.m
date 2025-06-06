function [y, theta] = evolutionary_model(pop, DE, X, Y)
% Author: Guilherme Gomes 
    m = DE.m; n = DE.n; k = DE.k; N = DE.N; 
    
% Bias considered    
if strcmp(DE.bias, 'Yes')
    
    % create a matrix to store exponent values
    mat = zeros(m, k, N);
    % create a matrix to store parameters values
    theta = zeros(m+1,N);
    % create a matrix to store model output for each generation
    y = zeros(numel(Y),N);
    % create a matrix to store Moore–Penrose pseudo-inverse matrix
    % (Giustolisi and Savic, 2006 - eq 16)
    P = zeros(m+1,m+1,N);
    
    for i = 1:N
        % reshape pouplation matrix (create N ES matrices)
       % mat(:,:,i) = (reshape(pop(i,:), m, k));
        mat(:,:,i) = (reshape(pop(i,:), k, m))';
        % create auxiliary variable for exponential operation
        R = zeros(numel(Y),k,m);
        % evaluate exponential operation for each parameter m
        for j = 1:m
            for kk = 1:k
                R(:,kk,j) = X(:,kk).^mat(j,kk,i);
            end
        end
        % create Z matrix
        Z(:,1:m) = prod(R,2);
        % add unitary column vector
        Z = [ones(size(Z,1),1) Z];
        % solution for parameters (least squares)
        theta(:,i) = pinv(Z)*Y;
        % Moore–Penrose pseudo-inverse matrix
        P(:,:,i) = pinv(Z)*pinv(Z)';

        % Now evaluate fitness of equations in the population
        % create auxiliary variable for sum operation
        M_sum = zeros(numel(Y),m);
        for j = 1:m
            M_sum(:,j) = theta(j+1,i) * Z(:,j+1);
        end
        % Remove all values in Z matrix, because it's not necessary
        Z = [];
        
        % model output for population i
        y(:,i) = theta(1,i) + sum(M_sum,2);   
    end

% Bias NOT considered    
else
        % create a matrix to store exponent values
    mat = zeros(m, k, N);
    % create a matrix to store parameters values
    theta = zeros(m,N);%zeros(m+1,N);
    % create a matrix to store model output for each generation
    y = zeros(numel(Y),N);
    % create a matrix to store Moore–Penrose pseudo-inverse matrix
    % (Giustolisi and Savic, 2006 - eq 16)
    P = zeros(m,m,N);%(m+1,m+1,N);
    
    for i = 1:N
        % reshape pouplation matrix (create N ES matrices)
       % mat(:,:,i) = (reshape(pop(i,:), m, k));
        mat(:,:,i) = (reshape(pop(i,:), k, m))';
        % create auxiliary variable for exponential operation
        R = zeros(numel(Y),k,m);
        % evaluate exponential operation for each parameter m
        for j = 1:m
            for kk = 1:k
                R(:,kk,j) = X(:,kk).^mat(j,kk,i);
            end
        end
        % create Z matrix
        Z(:,1:m) = prod(R,2);
        % add unitary column vector
        %Z = [ones(size(Z,1),1) Z];
        % solution for parameters (least squares)
        theta(:,i) = pinv(Z)*Y;
        % Moore–Penrose pseudo-inverse matrix
        P(:,:,i) = pinv(Z)*pinv(Z)';

        % Now evaluate fitness of equations in the population
        % create auxiliary variable for sum operation
        M_sum = zeros(numel(Y),m);
        for j = 1:m
            M_sum(:,j) = theta(j,i) * Z(:,j);
            %M_sum(:,j) = theta(j+1,i) * Z(:,j+1);
        end
        % Remove all values in Z matrix, because it's not necessary
        Z = [];
        
        % model output for population i
        y(:,i) = sum(M_sum,2);  
        %y(:,i) = theta(1,i) + sum(M_sum,2);  
    end
end