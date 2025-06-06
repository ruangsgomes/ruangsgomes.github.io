
% This function executes compromise programming multi-objective analisys
% Author: Guilherme Jose Cunha Gomes
% date: 04 / APR / 2020

function [Out_MO] = compromise (A, C, a)

% numero de criterios
n = numel(A(:,1));
% numero de alternativas
na = numel(A(1,:));
% checar se algum objetivo que possui os mesmos valores para todas as opcoes
B = []; 
for i = 1:n
    if round((sum(A(i,:))/numel(A(i,:))),4) == round(A(i,1),4) 
        % armazena o indice do criterio com mesmo valores
        B = [B; i];
    end
end
% volta à matriz de alternativas e criterios 
A(B,:) = [];
if isempty(A)
    %error('Options have the same values')
    Out_MO.Y = 1;  % stop execution
else

    %% daqui para baixo o programa é antigo
    % numero de criterios
    n = numel(A(:,1));
    % numero de alternativas
    na = numel(A(1,:));
    % estende a matriz para colocar os melhores valores e os 
    % piores valores
    c = nan(n,2);
    for i = 1:n
        if strcmp(C(i,:), 'min')
             c(i,1) = min(A(i,:));
             c(i,2) = max(A(i,:));
         else
             c(i,1) = max(A(i,:));
             c(i,2) = min(A(i,:));
        end
    end
    % inclui a matriz final
    A = [A c];

    % conjunto de valores de s
    s = [1 2 3];

    % matriz de proximidade final
    Ls = nan(na,numel(s));
    for i = 1:numel(s)
        for j = 1:na
            % computa somatorio para a alternativa
            % (loop nos criterios)
            som = nan(n,1);
            for k = 1:n
                som(k) =a(k)*(((A(k,end-1) - A(k,j))/(A(k,end-1) - A(k,end)))^s(i));
            end
    %         % verfica se ha valor NaN (divisao por zero e atribuir Inf)
    %         [row, ~] = find(isnan(som));
    %         som(row) = Inf;
            % armazena na matriz Ls
            Ls(j,i) = (sum(som))^(1/s(i));  
        end
    end

    % para s = infinito
    L_inf = nan(na,1);
    for j = 1:na
        % computa somatorio para a alternativa
        % (loop nos criterios)
        som = nan(n,1);
        for k = 1:n
            som(k) =a(k)*(((A(k,end-1) - A(k,j))/(A(k,end-1) - A(k,end))));
        end
            % armazena na matriz Ls
        L_inf(j) = max(som);  
    end

    % concanetar o resultado final
    L = [Ls L_inf];

    % melhor alternativa
    M = min(L(:,1:end));
    % determinacao da linha da melhor alternativa
    [row, column] = find(L==M);
    % a linha mais frequente significa a melhor alternativa
    Y = mode(row);
    %% output
    Out_MO.Y = Y;   % melhor alternativa
    Out_MO.L = L;   % matriz com valores de s
    Out_MO.A = A;  % nova matriz de alternativas
end