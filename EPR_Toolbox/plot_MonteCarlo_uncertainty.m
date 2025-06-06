function[p]= plot_MonteCarlo_uncertainty( n_gen, n_MonteCarlo, EPR_out_MonteCarlo_MODEGA, EPR_eval_MonteCarlo_MODEGA, BEST, Y_cal, Y_eval, Names )
%function[p]= plot_MonteCarlo_uncertainty(m_optimum, n_gen, n_MonteCarlo, EPR_out_MODEGA, EPR_out_MonteCarlo_MODEGA, EPR_eval_MonteCarlo_MODEGA )


% m_optimum                  = optimum number of terms
% n_gen                      = number of generations
% n_MonteCarlo               = number of Monte Carlo simulations
% EPR_out_MODEGA             = EPR estructure of a single optimization with m_optimum
% EPR_out_MonteCarlo_MODEGA  = EPR estructure of multiple optimization runs with m_optimum
% EPR_eval_MonteCarlo_MODEGA = EPR evaluation statistic for multiple runs


fax=16;
sv=1;

% Generate auxiliary variables
    SSE_mean = zeros(n_gen,1);
    DE_mean = zeros(n_gen,1);
    GA_mean = zeros(n_gen,1);
    R2_cal = zeros(n_MonteCarlo,1); R2_eval = zeros(n_MonteCarlo,1);
    RMSE_cal = zeros(n_MonteCarlo,1); RMSE_eval = zeros(n_MonteCarlo,1);   
    r_cal = zeros(n_MonteCarlo,1);r_eval = zeros(n_MonteCarlo,1);
    E_rel_cal = zeros(n_MonteCarlo,1); E_rel_eval = zeros(n_MonteCarlo,1);
    
    gen_R2_cal = zeros(n_MonteCarlo,1);gen_RMSE_cal = zeros(n_MonteCarlo,1);
    gen_MAE_cal = zeros(n_MonteCarlo,1);gen_E_rel_cal = zeros(n_MonteCarlo,1);
    gen_r_cal = zeros(n_MonteCarlo,1);
    
    gen_R2_eval = zeros(n_MonteCarlo,1);gen_RMSE_eval = zeros(n_MonteCarlo,1);
    gen_MAE_eval = zeros(n_MonteCarlo,1);gen_E_rel_eval = zeros(n_MonteCarlo,1);
    gen_r_eval = zeros(n_MonteCarlo,1);
    
% loop to fill auxiliary variables
for i =1:n_MonteCarlo
    R2_cal(i,1) = EPR_out_MonteCarlo_MODEGA(:,i).R2;
    R2_eval(i,1) = EPR_eval_MonteCarlo_MODEGA(:,i).R2;
    RMSE_cal(i,1) = EPR_out_MonteCarlo_MODEGA(:,i).RMSE;
    RMSE_eval(i,1) = EPR_eval_MonteCarlo_MODEGA(:,i).RMSE;
    r_cal(i,1) = EPR_out_MonteCarlo_MODEGA(:,i).r;
    r_eval(i,1) = EPR_eval_MonteCarlo_MODEGA(:,i).r; 
    E_rel_cal(i,1) = EPR_out_MonteCarlo_MODEGA(:,i).E_rel;
    E_rel_eval(i,1) = EPR_eval_MonteCarlo_MODEGA(:,i).E_rel; 
        for j =1:n_gen
           SSE_mean(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).SSE(j,1); % compute SSE    
           DE_mean(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).N_alg(j,1); % compute DE contribution from MODEGA
           GA_mean(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).N_alg(j,2); % compute GA contribution from MODEGA
           %gen_R2(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).generation_R2(j,1);
           %gen_RMSE(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).generation_RMSE(j,1);
           %gen_MAE(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).generation_MAE(j,1);
           %gen_E_rel(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).generation_E_rel(j,1);
           %gen_r(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).generation_r(j,1);
           
           % Grab statistics across generations for calibration
           gen_R2_cal(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).stats_cal(j,1);
           gen_RMSE_cal(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).stats_cal(j,2);
           gen_MAE_cal(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).stats_cal(j,3);
           gen_E_rel_cal(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).stats_cal(j,4);
           gen_r_cal(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).stats_cal(j,5);
           
           
           % Grab statistics across generations for evaluation
           gen_R2_eval(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).stats_eval(j,1);
           gen_RMSE_eval(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).stats_eval(j,2);
           gen_MAE_eval(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).stats_eval(j,3);
           gen_E_rel_eval(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).stats_eval(j,4);
           gen_r_eval(j,i) = EPR_out_MonteCarlo_MODEGA(:,i).stats_eval(j,5);
           
           
        end
end

% Compute 95% interval
[SSE.low, SSE.up] = confidence_intervals (SSE_mean');
[DE.low, DE.up] = confidence_intervals (DE_mean');
[GA.low, GA.up] = confidence_intervals (GA_mean');

%% Plot uncertainties
% parameters for figure and panel size
plotheight=6;%12
plotwidth=10;%8
subplotsx=1;
subplotsy=2;%3;   
leftedge=1.5;
rightedge=0.5;   
topedge=0.5;
bottomedge=1.5;
spacex=1.2;
spacey=0.7;
fontsize=12;    
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

% setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'units','normalized','outerposition',[0.1 0 .5 .9])

%--------------------- Plot A -----------------------------%
% define position and initial settings for the axis
%ax=axes('position',sub_pos{1,3},...
ax=axes('position',sub_pos{1,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
n_gen_vec=linspace(1,n_gen,n_gen);
P3 = Fill_Ranges(n_gen_vec,SSE.low,SSE.up, [0.5 0.5 0.5]); hold on
P4 = plot((1:sv:n_gen)',mean(SSE_mean,2),'r-','LineStyle','-','Linewidth',2);  
%P5 = plot((1:sv:n_gen)',EPR_out_MODEGA(:,m_optimum).SSE(1:sv:end),'k-','LineStyle','-','Linewidth',2);  

alpha(P3,.3);

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
%xlabel('Number of generations ','interpreter','latex','FontSize',fax);
ylabel('SSE','interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);
   
%set(gca,'YLim',[0 20000]);
set(gca,'xticklabel',{[]});
%leg1 = legend([P5, P4, P3],['MODEGA'],['Mean'],['95$\%$ uncertainty'],...
leg1 = legend([P4, P3],['Mean'],['95$\%$ uncertainty'],...
'location', 'northeast','Orientation','horizontal');
set(leg1,'Interpreter','latex');
text(0.9,0.65,'(a)','units','normalized','FontSize',fax,'interpreter','latex');   



%---------------------------- Plot b

%ax=axes('position',sub_pos{1,2},...
ax=axes('position',sub_pos{1,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
%n_gen=linspace(1,DE,DE);
P3 = Fill_Ranges(n_gen_vec,DE.low,DE.up, [0.0 0.0 0.5]); hold on
P4 = Fill_Ranges(n_gen_vec,GA.low,GA.up, [0.0 0.5 0.0]); hold on
P5 = plot((1:sv:n_gen)',mean(DE_mean,2),'b-','LineStyle','-','Linewidth',2);  
P6 = plot((1:sv:n_gen)',mean(GA_mean,2),'g-','LineStyle','-','Linewidth',2);  

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
xlabel('Number of generations ','interpreter','latex','FontSize',fax);
ylabel('$N_0$','interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);
alpha(P3,.3);
alpha(P4,.3);
set(gca,'YLim',[0 70]);
%set(gca,'xticklabel',{[]});
leg1 = legend([P3, P5, P4, P6],['DE 95$\%$ uncertainty'],['DE mean'],['GA 95$\%$ uncertainty'],['GA mean'],...
'location', 'north','Orientation','horizontal','NumColumns',2);
set(leg1,'Interpreter','latex');
text(0.9,0.75,'(b)','units','normalized','FontSize',fax,'interpreter','latex');   

%% BOX PLOT 

%markersizes for boxplot
msb = 5;

% boxcolor for eval dataset
bc_eval = [1 .1 .1];


% parameters for figure and panel size
plotheight=10;%10
plotwidth=10;%10
subplotsx=1;
subplotsy=4;%3;   
leftedge=1.5;
rightedge=0.5;   
topedge=0.5;
bottomedge=1.5;
spacex=1.2;
spacey=0.7;
fontsize=12;    
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

% setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'units','normalized','outerposition',[0.1 0 .5 .9])

%--------------------- Plot A -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,4},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
   
% Split the dataset between intervals to improve graphs 
% set intervals
n_points = n_gen/10; 

% For CALIBRATION data set
% Grab line 1 and 10 data
Gen_Y1 = gen_R2_cal(1:n_points-1:n_points,:);
% Grab from line 20 and every 10 lines posterior
Gen_Y2 = gen_R2_cal(n_points*2:n_points:end,:);
% Merge the dataset
Gen_Y_cal = [Gen_Y1; Gen_Y2];

% same for EVALUATION data set
Gen_Y1 = gen_R2_eval(1:n_points-1:n_points,:);
% Grab from line 20 and every 10 lines posterior
Gen_Y2 = gen_R2_eval(n_points*2:n_points:end,:);
% Merge the dataset
Gen_Y_eval = [Gen_Y1; Gen_Y2];

%Create vector to plot x axis
box_X = linspace(1,n_gen,n_gen); 
box_X=box_X';
% Grab line 1 and 10
Gen_X1 = box_X(1:n_points-1:n_points,:);
% Grab from line 20 every 10 lines
Gen_X2 = box_X(n_points*2:n_points:end,:);
% Merge datas
Gen_X = [Gen_X1 ; Gen_X2];
%boxchart(gen_R2'); hold on;


boxchart(Gen_Y_cal','MarkerSize',msb,'BoxFaceColor','k','MarkerColor','k'); hold on;
boxchart(Gen_Y_eval','MarkerSize',msb,'BoxFaceColor',bc_eval,'MarkerColor',bc_eval); hold on;

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
set(gca,'xticklabel',{[]});
%xlabel('Number of generations ','interpreter','latex','FontSize',fax);
ylabel('R$^2$','interpreter','latex','FontSize',fax);
grid on
%set(gca,'xticklabel',num2str(Gen_X));
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);
set(gca,'YLim',[0 1]);
%leg1 = legend([P5, P4, P3],['MODEGA'],['Mean'],['95$\%$ Uncertainty'],...
%'location', 'best','Orientation','horizontal');
%set(leg1,'Interpreter','latex');
%text(0.9,0.5,'(a)','units','normalized','FontSize',fax,'interpreter','latex');   

%--------------------- Plot B -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,3},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
   
% for CALIBRATION dataset        
% Split the dataset between intervals to improve graphs 
% Grab line 1 and 10 data
Gen_Y1 = gen_RMSE_cal(1:n_points-1:n_points,:);
% Grab from line 20 and every 10 lines posterior
Gen_Y2 = gen_RMSE_cal(n_points*2:n_points:end,:);
% Merge the dataset
Gen_Y_cal = [Gen_Y1; Gen_Y2];

% same for EVALUATION data set
Gen_Y1 = gen_RMSE_eval(1:n_points-1:n_points,:);
% Grab from line 20 and every 10 lines posterior
Gen_Y2 = gen_RMSE_eval(n_points*2:n_points:end,:);
% Merge the dataset
Gen_Y_eval = [Gen_Y1; Gen_Y2];

boxchart(Gen_Y_cal','MarkerSize',msb,'BoxFaceColor','k','MarkerColor','k'); hold on;
boxchart(Gen_Y_eval','MarkerSize',msb,'BoxFaceColor',bc_eval,'MarkerColor',bc_eval); hold on;

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
set(gca,'xticklabel',{[]});
%xlabel('Number of generations ','interpreter','latex','FontSize',fax);
ylabel('RMSE','interpreter','latex','FontSize',fax);
grid on
%set(gca,'xticklabel',num2str(Gen_X));
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);
set(gca,'YLim',[0 1]);
%leg1 = legend([P5, P4, P3],['MODEGA'],['Mean'],['95$\%$ Uncertainty'],...
%'location', 'best','Orientation','horizontal');
%set(leg1,'Interpreter','latex');
%text(0.9,0.5,'(b)','units','normalized','FontSize',fax,'interpreter','latex'); 

%--------------------- Plot C -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
   
% For CALIBRATION dataset        
% Split the dataset between intervals to improve graphs 
% Grab line 1 and 10 data
Gen_Y1 = gen_r_cal(1:n_points-1:n_points,:);
% Grab from line 20 and every 10 lines posterior
Gen_Y2 = gen_r_cal(n_points*2:n_points:end,:);
% Merge the dataset
Gen_Y_cal = [Gen_Y1; Gen_Y2];

% same for EVALUATION data set
Gen_Y1 = gen_r_eval(1:n_points-1:n_points,:);
% Grab from line 20 and every 10 lines posterior
Gen_Y2 = gen_r_eval(n_points*2:n_points:end,:);
% Merge the dataset
Gen_Y_eval = [Gen_Y1; Gen_Y2];

boxchart(Gen_Y_cal','MarkerSize',msb,'BoxFaceColor','k','MarkerColor','k'); hold on;
boxchart(Gen_Y_eval','MarkerSize',msb,'BoxFaceColor',bc_eval,'MarkerColor',bc_eval); hold on;

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
set(gca,'xticklabel',{[]});
%xlabel('Number of generations ','interpreter','latex','FontSize',fax);
ylabel('$r$','interpreter','latex','FontSize',fax);
grid on
%set(gca,'xticklabel',num2str(Gen_X));
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);
set(gca,'YLim',[0 1]);
%leg1 = legend([P5, P4, P3],['MODEGA'],['Mean'],['95$\%$ Uncertainty'],...
%'location', 'best','Orientation','horizontal');
%set(leg1,'Interpreter','latex');
%text(0.9,0.5,'(c)','units','normalized','FontSize',fax,'interpreter','latex'); 

%--------------------- Plot D -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
   
% For CALIBRATION dataset        
% Split the dataset between intervals to improve graphs 
% Grab line 1 and 10 data
Gen_Y1 = gen_E_rel_cal(1:n_points-1:n_points,:);
% Grab from line 20 and every 10 lines posterior
Gen_Y2 = gen_E_rel_cal(n_points*2:n_points:end,:);
% Merge the dataset
Gen_Y_cal = [Gen_Y1; Gen_Y2];


% same for EVALUATION data set
Gen_Y1 = gen_E_rel_eval(1:n_points-1:n_points,:);
% Grab from line 20 and every 10 lines posterior
Gen_Y2 = gen_E_rel_eval(n_points*2:n_points:end,:);
% Merge the dataset
Gen_Y_eval = [Gen_Y1; Gen_Y2];

P1 = boxchart(Gen_Y_cal','MarkerSize',msb,'BoxFaceColor','k','MarkerColor','k'); hold on;
P2 = boxchart(Gen_Y_eval','MarkerSize',msb,'BoxFaceColor',bc_eval,'MarkerColor',bc_eval); hold on;

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
xlabel('Number of generations ','interpreter','latex','FontSize',fax);
ylabel('$E_{rel}$','interpreter','latex','FontSize',fax);
grid on
set(gca,'xticklabel',num2str(Gen_X));
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);
set(gca,'YLim',[0 1]);
leg1 = legend([P1, P2],['Training'],['Testing'],...
'location', 'south','Orientation','horizontal');
set(leg1,'Interpreter','latex');
%text(0.9,0.5,'(d)','units','normalized','FontSize',fax,'interpreter','latex'); 




%% Box plot

% parameters for figure and panel size
%plotheight=10;
%plotwidth=10;
%subplotsx=2;
%subplotsy=2;   
%leftedge=1.5;
%rightedge=0.5;   
%topedge=0.5;
%bottomedge=1.5;
%spacex=1.2;
%spacey=1.5;
%fontsize=16;    
%sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

% setting the Matlab figure
%f=figure('visible','on');
%clf(f);
%set(gcf, 'PaperUnits', 'centimeters');
%set(gcf, 'PaperSize', [plotwidth plotheight]);
%set(gcf, 'PaperPositionMode', 'manual');
%set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
%set(gcf,'units','normalized','outerposition',[0.0 0.0 .7 .9])


plotheight=10;%10
plotwidth=5;%10
subplotsx=1;
subplotsy=4;%3;   
leftedge=1.5;
rightedge=0.5;   
topedge=0.5;
bottomedge=1.5;
spacex=1.2;
spacey=0.7;
fontsize=16;    
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

% setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'units','normalized','outerposition',[0.1 0 .2 .9])

%--------------------- Plot A -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,4},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');

sv = 1; % space vector

fax = fontsize;
color_vec = [1,1,1];
boxchart(R2_cal,'MarkerSize',msb,'BoxFaceColor','k','MarkerColor','k'); hold on
boxchart(R2_eval,'MarkerSize',msb,'BoxFaceColor',bc_eval,'MarkerColor',bc_eval); hold on

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );

set(gca,'xticklabel',{[]});
%ylabel('$R^2$','interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);
%set(gca,'xticklabel',num2str(n_gen));
   
% plot grid
grid on
% axes adjustment
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );

%text(0.05,0.90,'(a)','units','normalized','FontSize',fax,'interpreter','latex');  

%--------------------- Plot B -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,3},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
boxchart(RMSE_cal,'MarkerSize',msb,'BoxFaceColor','k','MarkerColor','k'); hold on
boxchart(RMSE_eval,'MarkerSize',msb,'BoxFaceColor',bc_eval,'MarkerColor',bc_eval); hold on

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );

%xlabel('Number of generations ','interpreter','latex','FontSize',fax);
%ylabel('\rm RMSE (\%)','interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{[]});
set(gca,'fontsize',fax);
   
% plot grid
grid on
% axes adjustment
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );

%text(0.05,0.90,'(b)','units','normalized','FontSize',fax,'interpreter','latex');  

%--------------------- Plot C -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
boxchart(r_cal,'MarkerSize',msb,'BoxFaceColor','k','MarkerColor','k'); hold on
boxchart(r_eval,'MarkerSize',msb,'BoxFaceColor',bc_eval,'MarkerColor',bc_eval); hold on

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );

%ylabel('$r$','interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);
set(gca,'xticklabel',{[]});
   
% plot grid
grid on
% axes adjustment
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );

%text(0.05,0.90,'(c)','units','normalized','FontSize',fax,'interpreter','latex'); 

%--------------------- Plot D -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
boxchart(E_rel_cal,'MarkerSize',msb,'BoxFaceColor','k','MarkerColor','k'); hold on
boxchart(E_rel_eval,'MarkerSize',msb,'BoxFaceColor',bc_eval,'MarkerColor',bc_eval); hold on

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );

%ylabel('$E_{rel}$','interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',num2str(n_gen));
set(gca,'fontsize',fax);
   
% plot grid
grid on
% axes adjustment
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );

%text(0.05,0.90,'(d)','units','normalized','FontSize',fax,'interpreter','latex'); 




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

% ref line
hline = refline(1,0);
set(hline,'Color','k','LineStyle','--');
xlabel(['Observed ' cellstr(Names(1))],'interpreter','latex','FontSize',fax);
ylabel(['Simulated ' cellstr(Names(1))],'interpreter','latex','FontSize',fax);
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

% ref line
hline = refline(1,0);
set(hline,'Color','k','LineStyle','--'); 

xlabel(['Observed ' cellstr(Names(1))],'interpreter','latex','FontSize',fax);
ylabel(['Simulated ' cellstr(Names(1))],'interpreter','latex','FontSize',fax);
grid on
set(gca,'fontsize',fax);
set(gca,'TickLabelInterpreter','latex')
text(0.05,0.90,'(b)','units','normalized','FontSize',fax,'interpreter','latex');   










p=fax;