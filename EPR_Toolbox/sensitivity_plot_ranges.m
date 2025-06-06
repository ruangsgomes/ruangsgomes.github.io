function [n]=plot_sensitivity_ranges(n,X,Y, X_cal, Y_cal, names)
% This function uses the number of Monte Carlo simulations (x)
% The sensitivity X matrix, and EPR 

% Fontsize
fax = 18;

% compute independent variables
n_x = size(X.mean,2);

% Number of parameters
for i=1:n_x
    figure;
    % Number of MonteCarlo runs
    for j = 1: n 
        X_value(:,j) = X.X(:,i,i);
        Y_value(:,j) = Y(j).Y_sensitivity(:,i);
        %plot(X_value(:,j), Y_value (:,j), 'k-','LineStyle','-','Linewidth',2); hold on
    end
    [Y_range.low, Y_range.up] = confidence_intervals (Y_value');
    %n_range=linspace(1,,n);
    Fill_Ranges(X.maxmin(:,i),Y_range.low,Y_range.up, [0.5 0.5 0.5]);hold on 
    P4 = plot(X.maxmin(:,i),mean(Y_value,2),'r--','LineStyle','-','Linewidth',1);  hold on
    P5 = plot(X_cal(:,i),Y_cal,'ko','MarkerSize',8);  hold on
    

    xlabel(cellstr(names(i+1)),'interpreter','latex','FontSize',fax);
    ylabel(cellstr(names(1)),'interpreter','latex','FontSize',fax);
    grid on
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'fontsize',fax);
    set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'LineWidth'   , 1        );

    if i==2
       set(gca,'XScale','log');
    end
    
    
end
    

   
%just for output
n.x=n_x;
n.y = n_y;