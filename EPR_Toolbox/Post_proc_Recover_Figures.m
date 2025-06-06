clear all; clc; close all;

load Run_EPR;


Names ={'$\sigma_1$' '$\beta^o$' '$\sigma_c$' 'R' '$c$' '$\phi$'};


%% PLOT

%---------------------  Plot uncertainty and average  --------------------%
[p]= plot_MonteCarlo_uncertainty(ndega, DE.last, n_MonteCarlo, EPR_out_MODEGA,EPR_out_MonteCarlo_MODEGA,EPR_eval_MonteCarlo_MODEGA);

%----------  Plot sensitivity uncertainty, mean and best model -----------%

% Plot sensitivity and output mean values of sensitivity
%[Y_mean_sensitivity] = plot_sensitivity_ranges(n_MonteCarlo,X_sensitivity, EPR_out_MonteCarlo_MODEGA, X_cal, Y_cal, Names);

% Function to select best model and plot 
[BEST]=model_selection(n_MonteCarlo,X_sensitivity,EPR_out_MonteCarlo_MODEGA, X_cal, Y_cal, Names);

%---------------------------- Best model graphs --------------------------%


% parameters for figure and panel size
plotheight=10;
plotwidth=10;
subplotsx=2;
subplotsy=2;   
leftedge=1.5;
rightedge=0.5;   
topedge=0.5;
bottomedge=1.5;
spacex=1.2;
spacey=1.5;
fontsize=14;    
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

% setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'units','normalized','outerposition',[0.1 0.1 .6 .9])

%--------------------- Plot A -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');

%% Post processing Compromise programming (R2)

% Now plot

fax = fontsize;
color_vec = [1,1,1];


P3 = plot(Y_cal,EPR_out_MonteCarlo_MODEGA(BEST).besty,'ks');
set(P3,'MarkerFaceColor','k',...
    'MarkerSize',8);

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );

set(gca,'YLim',[0 1000]);
set(gca,'XLim',[0 1000]);

% ref line
hline = refline(1,0);
set(hline,'Color','k');

xlabel('Observed $\sigma_1$','interpreter','latex','FontSize',fax);
ylabel('Simulated $\sigma_1$','interpreter','latex','FontSize',fax);

grid on
set(gca,'fontsize',fax);
set(gca,'TickLabelInterpreter','latex')

text(0.05,0.90,'(a)','units','normalized','FontSize',fax,'interpreter','latex');   

%--------------------- Plot B -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');


P5 = plot(Y_eval,EPR_eval_MonteCarlo_MODEGA(BEST).Y_eval,'ks');
set(P5,'MarkerFaceColor',color_vec,...
    'MarkerSize',8);
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );

set(gca,'YLim',[0 1000]);
set(gca,'XLim',[0 1000]);

% ref line
hline = refline(1,0);
set(hline,'Color','k'); 
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );


xlabel('Observed $\sigma_1$','interpreter','latex','FontSize',fax);
ylabel('Simulated $\sigma_1$','interpreter','latex','FontSize',fax);

grid on
set(gca,'fontsize',fax);
set(gca,'TickLabelInterpreter','latex')
text(0.05,0.90,'(b)','units','normalized','FontSize',fax,'interpreter','latex');   



%-------------------------------------------------------------------------%
%                             Sensitivity plot


% Number of parameters
for i=1:size(X_sensitivity.mean,2)
    figure;
    
    % Number of MonteCarlo runs
    for j = 1: n_MonteCarlo 
        Y_value(:,j) = EPR_out_MonteCarlo_MODEGA(j).Y_sensitivity(:,i);
    end
    
    % Get ranges
    [Y_range.low, Y_range.up] = confidence_intervals (Y_value');
    % Plot ranges
    P3 = Fill_Ranges(X_sensitivity.maxmin(:,i),Y_range.low,Y_range.up, [0.5 0.5 0.5]);hold on 
    alpha(P3,.3);
    % Plot mean values
    P4 = plot(X_sensitivity.maxmin(:,i),mean(Y_value,2),'r-','LineStyle','-','Linewidth',2);  hold on
    % Plot Best values
    P6 = plot(X_sensitivity.maxmin(:,i), EPR_out_MonteCarlo_MODEGA(BEST).Y_sensitivity(:,i),'k-','LineStyle','--','Linewidth',2);
   
    
    xlabel(cellstr(Names(i+1)),'interpreter','latex','FontSize',fax);
    ylabel(cellstr(Names(1)),'interpreter','latex','FontSize',fax);
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

    % Set log scale manually for a given variable
    if i==1
     %  leg1 = legend([ P3, P4, P6, P7, P8],['95$\%$ uncertainty'],['Mean'],...
%['MODEGA-SR1'], ['Middle model'],['Worst model'],'location', 'northeast');
       leg1 = legend([ P3, P4, P6],['95$\%$ uncertainty'],['Mean'],...
['Opt. Model'],'location', 'northeast');
       set(leg1,'Interpreter','latex');
    end   
    
   % if i==2
    %   set(gca,'XScale','log');
    %end
    set(gca,'YLim',[0 700]);
end
