function [foto] = Run_post_proc(X,DE,Output_DEGA)
% Run post processing

% Sigma 1 vs. Beta


%clear all; clc; close all;

foto=1;


fax = 14;


%Equation
DEGAMC_opt.m = DE.m;
DEGAMC_opt.k = DE.k;
DEGAMC_opt.bias = DE.bias;
Output_DEGAMC_opt.best_par = Output_DEGA.best_par;
Output_DEGAMC_opt.bestpop = Output_DEGA.bestpop;

% eq 4
%DEGAMC_opt.m = 3;
%DEGAMC_opt.k = 5;
%Output_DEGAMC_opt.best_par = [183.700875, -2.987641, 0.000027, 0.000603];
%Output_DEGAMC_opt.bestpop = [1.1,0,0,0,0 ; 0, 0.7, -0.1, 1.1, 2.6; 3, 0,0,0,0 ];

% eq 5
%DEGAMC_opt.m = 4;
%DEGAMC_opt.k = 5;
%Output_DEGAMC_opt.best_par = [160.213542, 0.000074, 121.824585, 0.322790, -1.237739];
%Output_DEGAMC_opt.bestpop = [0,0.7,0,1,2.4; 0,0,-3,0,0;2,0,0,0,0;1.7,0,0,0,0];

% eq 6
%DEGAMC_opt.m = 5;
%DEGAMC_opt.k = 5;
%Output_DEGAMC_opt.best_par = [159.041626,1.205461,0.000009, -2.974206, 85.799756, 0.000573];
%Output_DEGAMC_opt.bestpop = [0,1,0,0,0;0,0.5,0,1.1,3;1.1,0,0,0,0;0,0,-3,0,0;3,0,0,0,0];

%-----------------------HOEK and brown

%DEGAMC_opt.m = 4;
%DEGAMC_opt.k = 5;
%Output_DEGAMC_opt.best_par = [0.704042 0.002015 0.019532 3.965871 -0.030078];
%Output_DEGAMC_opt.bestpop = [2,1,0,0,0;0,0,0,1.70000000000000,0;0,0.900000000000002,-1,0,0.500000000000000;1.40000000000000,1,0,0,0];



%load DataBase_ConjuntoB_all.mat;
%X = DataBase_CojuntoB_all;
%parameters for figure and panel size
plotheight=10;
plotwidth=10;
subplotsx=5;
subplotsy=5;   
leftedge=1.5;
rightedge=0.5;   
topedge=0.5;
bottomedge=1.5;
spacex=0.6;%1.2;
spacey=1.;%1.5;
fontsize=16;    
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

% setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'units','normalized','outerposition',[0.0 0.0 .7 .9])

%--------------------- Plot 1 (gneiss Horino)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,5},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
%load Slate.mat;
%load Slate_ConjA_HB.mat;
%Slate = Slate_ConjA_HB;
%Slate = [Slate(:,1), Slate(:,2), Slate(:,3), Slate(:,4), Slate(:,5), Slate(:,6)];

X_GneissH = [X(1:60,1), X(1:60,2), X(1:60,3), X(1:60,4), X(1:60,5), X(1:60,6)];
X1 = linspace(0.1, 90,100)';
% X1 = deg2rad(X1);
% X1 = cos(4.*X1);

%figure;
for j = 1:6
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_GneissH(j,3), X_GneissH(j,4), X_GneissH(j,5), X_GneissH(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_GneissH(:,2),X_GneissH(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 800])
%xlim([-1 91])
ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Gneiss Horino','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );

%--------------------- Plot 2 (gneiss A Saroglou)-----------------------------%
% define position and initial settings for the axis
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,4},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
%load Slate.mat;
%load Slate_ConjA_HB.mat;
%Slate = Slate_ConjA_HB;
%Slate = [Slate(:,1), Slate(:,2), Slate(:,3), Slate(:,4), Slate(:,5), Slate(:,6)];

X_GneissA = [X(61:72,1), X(61:72,2), X(61:72,3), X(61:72,4), X(61:72,5), X(61:72,6)];


%figure;
for j = 1:4
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_GneissA(j,3), X_GneissA(j,4), X_GneissA(j,5), X_GneissA(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_GneissA(:,2),X_GneissA(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 400])
%xlim([-1 91])
ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Gneiss A','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 3 (Green River shale II)-----------------------------%
% define position and initial settings for the axis
% define position and initial settings for the axis
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,3},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_GRshaleII = [X(73:100,1), X(73:100,2), X(73:100,3), X(73:100,4), X(73:100,5), X(73:100,6)];


%figure;
for j = 1:4
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_GRshaleII(j,3), X_GRshaleII(j,4), X_GRshaleII(j,5), X_GRshaleII(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_GRshaleII(:,2),X_GRshaleII(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 600])
%xlim([-1 91])
ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Green River Shale II','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 4 (Moretown_phyllite)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_MTPhyllite = [X(101:121,1), X(101:121,2), X(101:121,3), X(101:121,4), X(101:121,5), X(101:121,6)];


%figure;
for j = 1:3
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_MTPhyllite(j,3), X_MTPhyllite(j,4), X_MTPhyllite(j,5), X_MTPhyllite(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_MTPhyllite(:,2),X_MTPhyllite(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 1000])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Moretown Phyllite','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );

%--------------------- Plot 5 (slate hoek)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_Slate = [X(122:161,1), X(122:161,2), X(122:161,3), X(122:161,4), X(122:161,5), X(122:161,6)];


%figure;
for j = 1:4
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_Slate(j,3), X_Slate(j,4), X_Slate(j,5), X_Slate(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_Slate(:,2),X_Slate(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 1000])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Slate Hoek','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );

%--------------------- Plot 6 (Penrhyn slate)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,5},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
X_PSlate = [X(162:203,1), X(162:203,2), X(162:203,3), X(162:203,4), X(162:203,5), X(162:203,6)];


%figure;
for j = 1:6
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_PSlate(j,3), X_PSlate(j,4), X_PSlate(j,5), X_PSlate(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_PSlate(:,2),X_PSlate(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 500])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Penhyrn Slate','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 7 (Martingsburg slate)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,4},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_MSlate = [X(204:238,1), X(204:238,2), X(204:238,3), X(204:238,4), X(204:238,5), X(204:238,6)];


%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_MSlate(j,3), X_MSlate(j,4), X_MSlate(j,5), X_MSlate(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_MSlate(:,2),X_MSlate(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 1000])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'MartingsBurg Slate','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 8 (Austin Slate)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,3},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_ASlate = [X(239:283,1), X(239:283,2), X(239:283,3), X(239:283,4), X(239:283,5), X(239:283,6)];


%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_ASlate(j,3), X_ASlate(j,4), X_ASlate(j,5), X_ASlate(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_ASlate(:,2),X_ASlate(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 1000])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Austin Slate','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );




%--------------------- Plot 9 (orthoquartzite dry Kumar 2006)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_OQuartz = [X(284:318,1), X(284:318,2), X(284:318,3), X(284:318,4), X(284:318,5), X(284:318,6)];


%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_OQuartz(j,3), X_OQuartz(j,4), X_OQuartz(j,5), X_OQuartz(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_OQuartz(:,2),X_OQuartz(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 1000])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Orthoquartzite dry Kumar 2006','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );



%--------------------- Plot 10 (phyllites dry Kumar 2006)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_Phyllite= [X(319:353,1), X(319:353,2), X(319:353,3), X(319:353,4), X(319:353,5), X(319:353,6)];


%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_Phyllite(j,3), X_Phyllite(j,4), X_Phyllite(j,5), X_Phyllite(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_Phyllite(:,2),X_Phyllite(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 500])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Phyllite dry Kumar 2006','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );



%--------------------- Plot 11 (slate dry Kumar 2006)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{3,5},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_Slate2006= [X(354:388,1), X(354:388,2), X(354:388,3), X(354:388,4), X(354:388,5), X(354:388,6)];


%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_Slate2006(j,3), X_Slate2006(j,4), X_Slate2006(j,5), X_Slate2006(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_Slate2006(:,2),X_Slate2006(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 500])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Slate dry Kumar 2006','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 12 (Tournemire shale)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{3,4},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
%load Slate.mat;
%load Slate_ConjA_HB.mat;
%Slate = Slate_ConjA_HB;
%Slate = [Slate(:,1), Slate(:,2), Slate(:,3), Slate(:,4), Slate(:,5), Slate(:,6)];

X_TourShale = [X(389:412,1), X(389:412,2), X(389:412,3), X(389:412,4), X(389:412,5), X(389:412,6)];

%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_TourShale(j,3), X_TourShale(j,4), X_TourShale(j,5), X_TourShale(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_TourShale(:,2),X_TourShale(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 800])
%xlim([-1 91])
ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Tournemire shale','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );

%--------------------- Plot 13 (limestone)-----------------------------%
% define position and initial settings for the axis
% define position and initial settings for the axis
ax=axes('position',sub_pos{3,3},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
%load Slate.mat;
%load Slate_ConjA_HB.mat;
%Slate = Slate_ConjA_HB;
%Slate = [Slate(:,1), Slate(:,2), Slate(:,3), Slate(:,4), Slate(:,5), Slate(:,6)];

X_lime = [X(413:447,1), X(413:447,2), X(413:447,3), X(413:447,4), X(413:447,5), X(413:447,6)];


%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_lime(j,3), X_lime(j,4), X_lime(j,5), X_lime(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_lime(:,2),X_lime(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 400])
%xlim([-1 91])
ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'limestone','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 14 (Green River shale I)-----------------------------%
% define position and initial settings for the axis
% define position and initial settings for the axis
% define position and initial settings for the axis
ax=axes('position',sub_pos{3,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_GRshaleI = [X(448:482,1), X(448:482,2), X(448:482,3), X(448:482,4), X(448:482,5), X(448:482,6)];


%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_GRshaleI(j,3), X_GRshaleI(j,4), X_GRshaleI(j,5), X_GRshaleI(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_GRshaleI(:,2),X_GRshaleI(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 600])
%xlim([-1 91])
ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Green River Shale I','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 15 (Carbonaceous Phyllite)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{3,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_CarbPhyllite = [X(483:517,1), X(483:517,2), X(483:517,3), X(483:517,4), X(483:517,5), X(483:517,6)];


%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_CarbPhyllite(j,3), X_CarbPhyllite(j,4), X_CarbPhyllite(j,5), X_CarbPhyllite(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_CarbPhyllite(:,2),X_CarbPhyllite(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 1000])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Carbonaceous Phyllite','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );

%--------------------- Plot 16 (Micaceous phyllite)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{4,5},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_MicPhyllite = [X(518:552,1), X(518:552,2), X(518:552,3), X(518:552,4), X(518:552,5), X(518:552,6)];


%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_MicPhyllite(j,3), X_MicPhyllite(j,4), X_MicPhyllite(j,5), X_MicPhyllite(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_MicPhyllite(:,2),X_MicPhyllite(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 1000])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Micaceous phyllite','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );

%--------------------- Plot 17 (Quartzite phyllite)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{4,4},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
X_QuartzPhyllite = [X(553:587,1), X(553:587,2), X(553:587,3), X(553:587,4), X(553:587,5), X(553:587,6)];


%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_QuartzPhyllite(j,3), X_QuartzPhyllite(j,4), X_QuartzPhyllite(j,5), X_QuartzPhyllite(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_QuartzPhyllite(:,2),X_QuartzPhyllite(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 500])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Quartzite phyllite','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 18 (Biotite schist)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{4,3},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_BioSchist = [X(588:617,1), X(588:617,2), X(588:617,3), X(588:617,4), X(588:617,5), X(588:617,6)];


%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_BioSchist(j,3), X_BioSchist(j,4), X_BioSchist(j,5), X_BioSchist(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_BioSchist(:,2),X_BioSchist(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 1000])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Biotite schist','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 19 (Chlorite schist)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{4,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_ChloSchist = [X(618:652,1), X(618:652,2), X(618:652,3), X(618:652,4), X(618:652,5), X(618:652,6)];


%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_ChloSchist(j,3), X_ChloSchist(j,4), X_ChloSchist(j,5), X_ChloSchist(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_ChloSchist(:,2),X_ChloSchist(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 1000])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Chlorite schist','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 20 (Quartzite schist)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{4,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_QuartSchist = [X(653:682,1), X(653:682,2), X(653:682,3), X(653:682,4), X(653:682,5), X(653:682,6)];

%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_QuartSchist(j,3), X_QuartSchist(j,4), X_QuartSchist(j,5), X_QuartSchist(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_QuartSchist(:,2),X_QuartSchist(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 1000])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Quartzite schist','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 21 (Quartz mica schist)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{5,5},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
X_QuartMicaSchist= [X(683:717,1), X(683:717,2), X(683:717,3), X(683:717,4), X(683:717,5), X(683:717,6)];

%figure;
for j = 1:5
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_QuartMicaSchist(j,3), X_QuartMicaSchist(j,4), X_QuartMicaSchist(j,5), X_QuartMicaSchist(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_QuartMicaSchist(:,2),X_QuartMicaSchist(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 500])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Quartz mica schist','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 22 (Layered rock)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{5,4},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
X_LayRock= [X(718:end,1), X(718:end,2), X(718:end,3), X(718:end,4), X(718:end,5), X(718:end,6)];

%figure;
for j = 1:6
    for i=1:size(X1,1)
      X_aux(i,:) = [X1(i,1), X_LayRock(j,3), X_LayRock(j,4), X_LayRock(j,5), X_LayRock(j,6)];
    end
    [s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);
    plot(X1(:,1),s_DEGA(:,1),'r-','MarkerFaceColor','k');hold on;
end

plot(X_LayRock(:,2),X_LayRock(:,1),'ko','MarkerFaceColor','k');hold on;

set(gca,'fontsize',fax)
ylim([0 500])
%xlim([-1 91])
%ylabel('$\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('$\beta$ (degree)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Layered rock','FontSize',fax,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
%set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

plotheight=10;
plotwidth=10;
subplotsx=5;
subplotsy=5;   
leftedge=1.5;
rightedge=0.5;   
topedge=0.5;
bottomedge=1.5;
spacex=0.6;%1.2;
spacey=1.;%1.5;
fontsize=16;    
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

% setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'units','normalized','outerposition',[0.0 0.0 .7 .9])

%--------------------- Plot 1 (gneiss Horino)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,5},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
%load Slate.mat;
%load Slate_ConjA_HB.mat;
%Slate = Slate_ConjA_HB;
%Slate = [Slate(:,1), Slate(:,2), Slate(:,3), Slate(:,4), Slate(:,5), Slate(:,6)];

X_GneissH = [X(1:60,1), X(1:60,2), X(1:60,3), X(1:60,4), X(1:60,5), X(1:60,6)];

X_aux = X_GneissH(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_GneissH(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.GneissH]=statistic_model(X_GneissH(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.GneissH.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.GneissH.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 800])
xlim([0 800])
ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Gneiss Horino','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );

%--------------------- Plot 2 (gneiss A Saroglou) -----------------------------%
% define position and initial settings for the axis
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,4},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
%load Slate.mat;
%load Slate_ConjA_HB.mat;
%Slate = Slate_ConjA_HB;
%Slate = [Slate(:,1), Slate(:,2), Slate(:,3), Slate(:,4), Slate(:,5), Slate(:,6)];

X_GneissA = [X(61:72,1), X(61:72,2), X(61:72,3), X(61:72,4), X(61:72,5), X(61:72,6)];

X_aux = X_GneissA(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_GneissA(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.GneissA]=statistic_model(X_GneissA(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.GneissA.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.GneissA.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 500])
xlim([0 500])
ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Gneiss A Saroglou','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 3 (Green River shale II) -----------------------------%
% define position and initial settings for the axis
% define position and initial settings for the axis
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,3},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_GRshaleII = [X(73:100,1), X(73:100,2), X(73:100,3), X(73:100,4), X(73:100,5), X(73:100,6)];

X_aux = X_GRshaleII(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_GRshaleII(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.GRshaleII]=statistic_model(X_GRshaleII(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.GRshaleII.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.GRshaleII.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 600])
xlim([0 600])
ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Green Rivel shale II','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 4 (Moretown_phyllite) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
X_MTPhyllite = [X(101:121,1), X(101:121,2), X(101:121,3), X(101:121,4), X(101:121,5), X(101:121,6)];

X_aux = X_MTPhyllite(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_MTPhyllite(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.MTPhyllite]=statistic_model(X_MTPhyllite(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.MTPhyllite.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.MTPhyllite.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 1000])
xlim([0 1000])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Moretown Phyllite','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );




%--------------------- Plot 5 (slate hoek) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_Slate = [X(122:161,1), X(122:161,2), X(122:161,3), X(122:161,4), X(122:161,5), X(122:161,6)];


X_aux = X_Slate(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_Slate(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.Slate]=statistic_model(X_Slate(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.Slate.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.Slate.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 1000])
xlim([0 1000])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Slate Hoek','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );



%--------------------- Plot 6 (Penrhyn slate) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,5},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_PSlate = [X(162:203,1), X(162:203,2), X(162:203,3), X(162:203,4), X(162:203,5), X(162:203,6)];


X_aux = X_PSlate(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_PSlate(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.PSlate]=statistic_model(X_PSlate(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.PSlate.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.PSlate.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 500])
xlim([0 500])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Penhyrn Slate','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );



%--------------------- Plot 7 (Martingsburg slate) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,4},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_MSlate = [X(204:238,1), X(204:238,2), X(204:238,3), X(204:238,4), X(204:238,5), X(204:238,6)];

X_aux = X_MSlate(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_MSlate(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.MSlate]=statistic_model(X_MSlate(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.MSlate.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.MSlate.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 800])
xlim([0 800])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Martingsburg slate','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );



%--------------------- Plot 8 (Austin Slate) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,3},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_ASlate = [X(239:283,1), X(239:283,2), X(239:283,3), X(239:283,4), X(239:283,5), X(239:283,6)];

X_aux = X_ASlate(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_ASlate(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.ASlate]=statistic_model(X_ASlate(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.ASlate.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.ASlate.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 1000])
xlim([0 1000])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Austin Slate','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );

%--------------------- Plot 9 (orthoquartzite dry Kumar 2006) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_OQuartz = [X(284:318,1), X(284:318,2), X(284:318,3), X(284:318,4), X(284:318,5), X(284:318,6)];

X_aux = X_OQuartz(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_OQuartz(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.OQuartz]=statistic_model(X_OQuartz(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.OQuartz.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.OQuartz.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')

set(gca,'fontsize',fax)
ylim([0 1000])
xlim([0 1000])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'orthoquartzite dry Kumar 2006','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );

%--------------------- Plot 10 (phyllites dry Kumar 2006) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_Phyllite= [X(319:353,1), X(319:353,2), X(319:353,3), X(319:353,4), X(319:353,5), X(319:353,6)];

X_aux = X_Phyllite(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_Phyllite(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.Phyllite]=statistic_model(X_Phyllite(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.Phyllite.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.Phyllite.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 500])
xlim([0 500])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'phyllites dry Kumar 2006','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );



%--------------------- Plot 11 (slate dry Kumar 2006) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{3,5},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_Slate2006= [X(354:388,1), X(354:388,2), X(354:388,3), X(354:388,4), X(354:388,5), X(354:388,6)];


X_aux = X_Slate2006(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_Slate2006(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.Slate2006]=statistic_model(X_Slate2006(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.Slate2006.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.Slate2006.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 500])
xlim([0 500])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'slate dry Kumar 2006','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );



%--------------------- Plot 12 (Tournemire shale)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{3,4},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
%load Slate.mat;
%load Slate_ConjA_HB.mat;
%Slate = Slate_ConjA_HB;
%Slate = [Slate(:,1), Slate(:,2), Slate(:,3), Slate(:,4), Slate(:,5), Slate(:,6)];

X_TourShale = [X(389:412,1), X(389:412,2), X(389:412,3), X(389:412,4), X(389:412,5), X(389:412,6)];

X_aux = X_TourShale(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_TourShale(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.TourShale]=statistic_model(X_TourShale(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.TourShale.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.TourShale.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 800])
xlim([0 800])
ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Tournemire shale','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );

%--------------------- Plot 13 (limestone) -----------------------------%
% define position and initial settings for the axis
% define position and initial settings for the axis
ax=axes('position',sub_pos{3,3},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
%load Slate.mat;
%load Slate_ConjA_HB.mat;
%Slate = Slate_ConjA_HB;
%Slate = [Slate(:,1), Slate(:,2), Slate(:,3), Slate(:,4), Slate(:,5), Slate(:,6)];

X_lime = [X(413:447,1), X(413:447,2), X(413:447,3), X(413:447,4), X(413:447,5), X(413:447,6)];

X_aux = X_lime(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_lime(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.lime]=statistic_model(X_lime(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.lime.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.lime.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 500])
xlim([0 500])
ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'limestone','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 14 (Green River shale I) -----------------------------%
% define position and initial settings for the axis
% define position and initial settings for the axis
% define position and initial settings for the axis
ax=axes('position',sub_pos{3,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_GRshaleI = [X(448:482,1), X(448:482,2), X(448:482,3), X(448:482,4), X(448:482,5), X(448:482,6)];

X_aux = X_GRshaleI(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_GRshaleI(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.GRshaleI]=statistic_model(X_GRshaleI(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.GRshaleI.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.GRshaleI.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 600])
xlim([0 600])
ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Green Rivel shale I','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 15 (Carbonaceous Phyllite) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{3,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
X_CarbPhyllite = [X(483:517,1), X(483:517,2), X(483:517,3), X(483:517,4), X(483:517,5), X(483:517,6)];

X_aux = X_CarbPhyllite(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_CarbPhyllite(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.CarbPhyllite]=statistic_model(X_CarbPhyllite(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.CarbPhyllite.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.CarbPhyllite.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 1000])
xlim([0 1000])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Carbonaceous Phyllite','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 16 (Micaceous phyllite) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{4,5},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
X_MicPhyllite = [X(518:552,1), X(518:552,2), X(518:552,3), X(518:552,4), X(518:552,5), X(518:552,6)];

X_aux = X_MicPhyllite(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_MicPhyllite(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.MicPhyllite]=statistic_model(X_MicPhyllite(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.MicPhyllite.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.MicPhyllite.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 1000])
xlim([0 1000])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Micaceous phyllite','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%--------------------- Plot 17 (Quartzite phyllite) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{4,4},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
X_QuartzPhyllite = [X(553:587,1), X(553:587,2), X(553:587,3), X(553:587,4), X(553:587,5), X(553:587,6)];

X_aux = X_QuartzPhyllite(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_QuartzPhyllite(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.QuartzPhyllite]=statistic_model(X_QuartzPhyllite(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.QuartzPhyllite.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.QuartzPhyllite.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 500])
xlim([0 500])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Quartzite phyllite','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );

%--------------------- Plot 18 (Biotite schist) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{4,3},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        
X_BioSchist = [X(588:617,1), X(588:617,2), X(588:617,3), X(588:617,4), X(588:617,5), X(588:617,6)];

X_aux = X_BioSchist(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_BioSchist(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.BioSchist]=statistic_model(X_BioSchist(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.BioSchist.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.BioSchist.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')

set(gca,'fontsize',fax)
ylim([0 800])
xlim([0 800])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Biotite schist','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );


%---------------------  Plot 19 (Chlorite schist) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{4,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_ChloSchist = [X(618:652,1), X(618:652,2), X(618:652,3), X(618:652,4), X(618:652,5), X(618:652,6)];

X_aux = X_ChloSchist(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_ChloSchist(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.ChloSchist]=statistic_model(X_ChloSchist(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.ChloSchist.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.ChloSchist.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 1000])
xlim([0 1000])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Chlorite schist','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );

%--------------------- Plot 20 (Quartzite schist) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{4,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_QuartSchist = [X(653:682,1), X(653:682,2), X(653:682,3), X(653:682,4), X(653:682,5), X(653:682,6)];

X_aux = X_QuartSchist(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_QuartSchist(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.QuartSchist]=statistic_model(X_QuartSchist(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.QuartSchist.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.QuartSchist.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')

set(gca,'fontsize',fax)
ylim([0 1000])
xlim([0 1000])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Quartzite schist','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );

%---------------------Plot 21 (Quartz mica schist) -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{5,5},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_QuartMicaSchist= [X(683:717,1), X(683:717,2), X(683:717,3), X(683:717,4), X(683:717,5), X(683:717,6)];

X_aux = X_QuartMicaSchist(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_QuartMicaSchist(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.QuartMicaSchist]=statistic_model(X_QuartMicaSchist(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.QuartMicaSchist.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.QuartMicaSchist.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 500])
xlim([0 500])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
%xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Quartz mica schist','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );



%--------------------- Plot 22 (Layered rock)-----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{5,4},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
        

X_LayRock= [X(718:end,1), X(718:end,2), X(718:end,3), X(718:end,4), X(718:end,5), X(718:end,6)];

X_aux = X_LayRock(:,2:end);

[s_DEGA] = evaluation_DE_regression (DEGAMC_opt, Output_DEGAMC_opt, X_aux);

plot(X_LayRock(:,1),s_DEGA(:,1),'ko','MarkerFaceColor','w');hold on;
plot(linspace(0,1500,1500),linspace(0,1500,1500),'k--');hold on;

[Stats.LayRock]=statistic_model(X_LayRock(:,1),s_DEGA(:,1));
text(0.63, .25,['R$^2$=',num2str(Stats.LayRock.R2,'%.2f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')
text(0.5, .1,['RMSE=',num2str(Stats.LayRock.RMSE,'%.1f')],'FontSize',fax-4,'Units','normalized','interpreter','latex')


set(gca,'fontsize',fax)
ylim([0 500])
xlim([0 500])
%ylabel('Pred. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
xlabel('Obs. $\sigma_{\rm 1}$ (MPa)','interpreter','latex','FontSize',fax);
text(0.05, .85,'Layered rock','FontSize',fax-2,'Units','normalized','interpreter','latex')
grid on

set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

set(gca, ...
  'Fontsize'    , fax , ...  
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'     , ...
  'LineWidth'   , 1         );