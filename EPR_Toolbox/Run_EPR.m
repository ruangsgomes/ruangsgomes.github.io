
% --- MultiObjective Differential Evolution Genetic Algorithm EPR ------- %
%                                                                         %
% MultiObjective Differential Evolution Genetic Algorithm EPR (EPR-MODEGA)%
% executes a dual search based evolutionary polynomial regression (EPR)   %
% with self-adaptive offspring creation and compromise programming model  % 
% selection                                                               %
%                                                                         %
%                                                                         %
% EPR-MODEGA developed by Guilherme Jose Cunha Gomes                      %
%                                                                         %
% ------------------------------------------------------------------------%
%                                                                         %
% The EPR-MODEGA code has been described in:                              %
%                                                                         %
%  Gomes, G.J.C; Gomes, R. Vargas Jr, E.A. (2020) A dual search based EPR %
%  with self-adaptive offspring creation and compromise programming model %
%  selection, Eng. Comp., XXXX.                                           %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
%     Copyright (C) 2019-2029  the author                                 %
%                                                                         %
%     This program is free software: you can modify it under the terms of %
%     the GNU General Public License as published by the Free Software    %
%     Foundation, either version 3 of the License, or (at your option)    %
%     later version.                                                      %
%                                                                         %
%     This program is distributed in the hope that it will be useful,     %
%     but WITHOUT ANY WARRANTY; without even the implied warranty of      %
%     without even the implied warranty of MERCHANTABILITY or FITNESS FOR %
%     A PARTICULAR PURPOSE. See the GNU General Public License for more   %
%     details.                                                            %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% MATLAB code written by Guilherme J. C. Gomes,                           %                             
%                                                                         %
%                        guilhermejcg@ufop.edu.br                         %
% Version 0.5: June 2020                                                  %                                                 
%                                                                         % 
%-------------------------------------------------------------------------%

% For more information please read: 

% Giustolisi and Savic (2006) A symbolic data-driven technique
% based on evlutionary polynomial regression. J. Hydroinformatics, 8(3),
% 207-222.
%
% A. Ahangar-Asr, A. Faramarzi, N. Mottaghifard, A.A. Javadi
% (2011) Modeling of permeability and compaction characteristics of soils 
% using evolutionary polynomial regression. Computers & Geosciences

% Clear memory and close matlab Figures
clc; clear all; close all;

% Clear all text files
delete output_EPR.txt
delete Model_Selection.txt
delete MonteCarlo_EPR.txt
delete SplitSampling.txt

%% 

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%                              EPR INPUTS
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% Insert dependent variable and independent variables names for plot
Names ={'K' '$\Phi$' 'F' 'dPT' 'dPore' 'dGrain'};

% Dataset will be provided already differentiated by train and test?
% ('Yes' or 'No')
Data_divided = 'Yes';

% If Data_divided = 'No' then select training contribution
if strcmp(Data_divided, 'No')
    % Set contribution in decimals (e.g. 75% = 0.75)
    P = 0.75;
end

% DEFAULT - Reads excel according to statement above
if  strcmp(Data_divided, 'Yes')
    % Read training/calibration dataset
    filename = 'Data_train.xlsx';
    Data1 = xlsread(filename); 
    Y_cal = Data1(2:end,1); X_cal = Data1(2:end,2:end);
    
    % Read testing/evaluation dataset
    filename = 'Data_test.xlsx';
    Data2 = xlsread(filename);
    Y_eval = Data2(2:end,1); X_eval = Data2(2:end,2:end); 
    
    % Compile entire dataset
    Y = [Y_cal; Y_eval]; X = [X_cal; X_eval]; 
    
    % Training Percentage
    P = round(size(Data1,1)/(size(Data1,1) + size(Data2,1)),2);
else
    % Read database
    filename = 'Data.xlsx';
    Data = xlsread(filename); Y = Data(2:end,1); X = Data (2:end,2:end);
    % Split dataset to 75% train and 25% test with similar mean value
    [Y_cal, X_cal, Y_eval, X_eval] = split_dataset(Y,X,P);
end

% Choose optional Bias ('Yes' or 'No');
DE.bias = 'Yes';
GA.bias = DE.bias;

% ----------------------- EVOLUTIONARY ALGORITHMS ----------------------- %
% Differential Evolution inputs
DE.last = 300;             % number of generations
DE.cr = 0.7;               % crossover probability rate
DE.gamma = 1;              % jump rate gamma
DE.delta = 1;              % number of elements selected in the algorithm      
DE.exponents = [-2:0.1:2]; % exponents
DE.max = 5;                % maximum number of EPR terms (m-value)
DE.k = size(X_cal,2);      % number of independent variables
DE.n = size(Y_cal,1);      % number of target values
DE.pop = 20;               % number of population multiplied by poly length

% Genetic Algorithm inputs
GA.last = DE.last;         % number of generations
GA.cpr = 0.4;              % crossover probability rate
GA.mutrate = 0.1;          % mutation rate
GA.k = size(X_cal,2);      % number of independent variables
GA.n = size(Y_cal,1);      % number of target values
GA.exponents = [-2:0.1:2]; % exponents   
    
%------------------------ COMPROMISE PROGRAMMING -------------------------%
% define statistic metrics to evaluate in the compromise-based
% multiobjective analysis
C = ['max'; % R2 (coeffient of determination)
     'min'; % RMSE (root-mean-square error)
     'max'; % r (coeffient of correlation)
     'max'; % E_rel
     'min'; % Num parametros
    ];
a = ones(size(C,1),1);      % objective wheights   

% create matrix to store data
CC =  nan(size(a,1),DE.max); %  MODE
DD =  nan(size(a,1),DE.max); %  MOGA
EE =  nan(size(a,1),DE.max); %  MODEGA

% ---------------------- SENSITIVITY-DRIVEN ----------------------------- %
% Chose if MODEGA-SD will be performed ('Yes' or 'No');
MODEGASD = 'Yes';
% Choose number of MonteCarlo
n_MonteCarlo = 10;

% ------------------------ CROSS-VALIDATION ----------------------------- %
% Choose if Cross-validated MODEGA-SD will be performed ('Yes' or 'No');
CrossVal = 'Yes';
% if 'Yes' then enter number of split samples for Cross-validation
n_split = 5;

% DEFAULT - If Cross validation 'No' stop code after 1 Monte Carlo
if  strcmp(CrossVal, 'No')
    n_split = 1; % Do not change this!
end

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
%                              END OF  EPR INPUTS
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

%%

hw = waitbar(0,'Running Multi-objective EPR...');
tic % start timer

% EPR optimization process
    for j = 1:DE.max
        
        %update waitbar
        waitbar(j/DE.max);
        
       % MODE Algorithm
       DE.m = j;         % number of EPR terms (excluding bias)
       DE.N = j*DE.pop;      % Population size
       [Output_MODE] = MODE_algorithm (DE, Y_cal, X_cal);
       
       % save data in a structure variable
        EPR_out_MODE(:,j).SSE = Output_MODE.SEE;           % best SSE for each generation 
        EPR_out_MODE(:,j).besty = Output_MODE.besty;       % best simulation
        EPR_out_MODE(:,j).bestpop = Output_MODE.bestpop;   % best population
        EPR_out_MODE(:,j).best_par = Output_MODE.best_par; % best EPR parameters
        EPR_out_MODE(:,j).R2 = Output_MODE.R2;             % R2
        EPR_out_MODE(:,j).RMSE = Output_MODE.RMSE;         % RMSE
        EPR_out_MODE(:,j).MAE = Output_MODE.MAE;           % MAE
        EPR_out_MODE(:,j).E_rel = Output_MODE.E_rel;       % Relative efficiency
        EPR_out_MODE(:,j).r = Output_MODE.r;               % Pearson's correlation 
        EPR_out_MODE(:,j).OF = Output_MODE.OF;             % best OF-values for each generation (SSE and num. variables in the model)
        
        % run EPR equation with evaluation data set
        [Y_MODE_eval] = evaluation_DE_regression (DE, Output_MODE, X_eval);
        EPR_eval_MODE(:,j).Y_eval = Y_MODE_eval;

        % compute statistics MOGA
        [Out_stat] = statistic_model (Y_eval, Y_MODE_eval);
        EPR_eval_MODE(:,j).R2 = Out_stat.R2;
        EPR_eval_MODE(:,j).RMSE = Out_stat.RMSE;
        EPR_eval_MODE(:,j).MAE = Out_stat.MAE;
        EPR_eval_MODE(:,j).E_rel = Out_stat.E_rel;
        EPR_eval_MODE(:,j).r = Out_stat.r;
       save_data (Output_MODE,DE,j,'MODE');

       % MOGA Algorithm
       GA.m = j;         % number of EPR terms (excluding bias)
       GA.N = j*20;      % Population size
       [Output_MOGA] = MOGA_algorithm (GA, Y_cal, X_cal);
       
       % save data in a structure variable
        EPR_out_MOGA(:,j).SSE = Output_MOGA.SEE;           % best SSE for each generation 
        EPR_out_MOGA(:,j).besty = Output_MOGA.besty;       % best simulation
        EPR_out_MOGA(:,j).bestpop = Output_MOGA.bestpop;   % best population
        EPR_out_MOGA(:,j).best_par = Output_MOGA.best_par; % best EPR parameters
        EPR_out_MOGA(:,j).R2 = Output_MOGA.R2;             % R2
        EPR_out_MOGA(:,j).RMSE = Output_MOGA.RMSE;         % RMSE
        EPR_out_MOGA(:,j).MAE = Output_MOGA.MAE;           % MAE
        EPR_out_MOGA(:,j).E_rel = Output_MOGA.E_rel;       % Relative efficiency
        EPR_out_MOGA(:,j).r = Output_MOGA.r;               % Pearson's correlation 
        EPR_out_MOGA(:,j).OF = Output_MOGA.OF;             % best OF-values for each generation (SSE and num. variables in the model)
        
        % run EPR equation with evaluation data set
        [Y_MOGA_eval] = evaluation_DE_regression (GA, Output_MOGA, X_eval);
        EPR_eval_MOGA(:,j).Y_eval = Y_MOGA_eval;

        % compute statistics MOGA
        [Out_stat] = statistic_model (Y_eval, Y_MOGA_eval);
        EPR_eval_MOGA(:,j).R2 = Out_stat.R2;
        EPR_eval_MOGA(:,j).RMSE = Out_stat.RMSE;
        EPR_eval_MOGA(:,j).MAE = Out_stat.MAE;
        EPR_eval_MOGA(:,j).E_rel = Out_stat.E_rel;
        EPR_eval_MOGA(:,j).r = Out_stat.r;       
        save_data (Output_MOGA,GA,j,'MOGA');
      
       % MODEGA Algorithm
       DE.m = j;         % number of EPR terms (excluding bias)
       DE.N = j*DE.pop;      % Population size
       [Output_MODEGA] = MODEGA_algorithm (DE, GA, Y_cal, X_cal);
       % save data in a structure variable
        EPR_out_MODEGA(:,j).SSE = Output_MODEGA.SEE;           % best SSE for each generation 
        EPR_out_MODEGA(:,j).besty = Output_MODEGA.besty;       % best simulation
        EPR_out_MODEGA(:,j).bestpop = Output_MODEGA.bestpop;   % best population
        EPR_out_MODEGA(:,j).best_par = Output_MODEGA.best_par; % best EPR parameters
        EPR_out_MODEGA(:,j).R2 = Output_MODEGA.R2;             % R2
        EPR_out_MODEGA(:,j).RMSE = Output_MODEGA.RMSE;         % RMSE
        EPR_out_MODEGA(:,j).MAE = Output_MODEGA.MAE;           % MAE
        EPR_out_MODEGA(:,j).E_rel = Output_MODEGA.E_rel;       % Relative efficiency
        EPR_out_MODEGA(:,j).r = Output_MODEGA.r;               % Pearson's correlation 
        EPR_out_MODEGA(:,j).OF = Output_MODEGA.OF;             % best OF-values for each generation (SSE and num. variables in the model)
        EPR_out_MODEGA(:,j).N_alg = Output_MODEGA.N_alg;             % number of offspring
        % run EPR equation with evaluation data set
        [Y_MODEGA_eval] = evaluation_DE_regression (DE, Output_MODEGA, X_eval);
        EPR_eval_MODEGA(:,j).Y_eval = Y_MODEGA_eval;

        % compute statistics MOGA
        [Out_stat] = statistic_model (Y_eval, Y_MODEGA_eval);
        EPR_eval_MODEGA(:,j).R2 = Out_stat.R2;
        EPR_eval_MODEGA(:,j).RMSE = Out_stat.RMSE;
        EPR_eval_MODEGA(:,j).MAE = Out_stat.MAE;
        EPR_eval_MODEGA(:,j).E_rel = Out_stat.E_rel;
        EPR_eval_MODEGA(:,j).r = Out_stat.r;       
        % save data in a structure variable
       save_data (Output_MODEGA,DE,j,'MODEGA');
     
       % MODE Algorithm
       CC(:,j)  = [EPR_out_MODE(:,j).R2;                % R2
                EPR_out_MODE(:,j).RMSE;                 % RMSE
                EPR_out_MODE(:,j).r;                    % r
                EPR_out_MODE(:,j).E_rel;                % E_rel
                numel(EPR_out_MODE(:,j).best_par) ];    % Num parameters

       % MOGA Algorithm
        DD(:,j)  = [EPR_out_MOGA(:,j).R2;                % R2
                EPR_out_MOGA(:,j).RMSE;                 % RMSE
                EPR_out_MOGA(:,j).r;                    % r
                EPR_out_MOGA(:,j).E_rel;                % E_rel
                numel(EPR_out_MOGA(:,j).best_par) ];    % Num parameters
      
       % MODEGA Algorithm
        EE(:,j)  = [EPR_out_MODEGA(:,j).R2;                % R2
                EPR_out_MODEGA(:,j).RMSE;                 % RMSE
                EPR_out_MODEGA(:,j).r;                    % r
                EPR_out_MODEGA(:,j).E_rel;                % E_rel
                numel(EPR_out_MODEGA(:,j).best_par) ];    % Num parameters
            
    end
     

% run Compromise MODE programming
[Out_MO] = compromise (CC, C, a);
% optimum number of EPR terms (nde) and Compromise matrix (L_de)
nde_MO = Out_MO.Y; L_de_MO = Out_MO.L;
% run Compromise MOGA programming
[Out_MO] = compromise (DD, C, a);
% optimum number of EPR terms (nde) and Compromise matrix (L_de)
nga_MO = Out_MO.Y; L_ga_MO = Out_MO.L;
% run Compromise DEGA programming
[Out_MO] = compromise (EE, C, a);
% optimum number of EPR terms (nde) and Compromise matrix (L_de)
ndega = Out_MO.Y; L_dega = Out_MO.L;

% run Compromise programming
 A = [CC DD EE ];    % concanetate all models

%% Tables with results    
% output file
fileID = fopen('output_EPR.txt','a');
% Create Table with dataset
Tab = [min(Y)' max(Y)' mean(Y)' std(Y)'; min(X)' max(X)' mean(X)' std(X)'];
fprintf(fileID,'%12s\r\n',' ');
fprintf(fileID,'%12s\r\n','---------- Tabelas publicacao ----------');
fprintf(fileID,'%12s\r\n','Statistics of variables used in the database:');
fprintf(fileID,'%12s\r\n','min, max, mean and standard deviation');
fprintf(fileID,'%6.4f %6.4f %6.4f %6.4f\n',Tab');
fprintf(fileID,'%12s\r\n',' ');
fprintf(fileID,'%12s\r\n','Summary of five indicators for the optimal model on training and testing data:');


% Table with calibration and validation data
Tab_DE = [ndega+1 ndega+1 nga_MO+1 nga_MO+1 nde_MO+1 nde_MO+1;...
    EPR_out_MODEGA(ndega).R2 EPR_eval_MODEGA(ndega).R2 EPR_out_MOGA(nga_MO).R2 EPR_eval_MOGA(nga_MO).R2 EPR_out_MODE(nde_MO).R2 EPR_eval_MODE(nde_MO).R2;...
    EPR_out_MODEGA(ndega).RMSE EPR_eval_MODEGA(ndega).RMSE EPR_out_MOGA(nga_MO).RMSE EPR_eval_MOGA(nga_MO).RMSE EPR_out_MODE(nde_MO).RMSE EPR_eval_MODE(nde_MO).RMSE;...
    EPR_out_MODEGA(ndega).r EPR_eval_MODEGA(ndega).r EPR_out_MOGA(nga_MO).r EPR_eval_MOGA(nga_MO).r EPR_out_MODE(nde_MO).r EPR_eval_MODE(nde_MO).r;...
    EPR_out_MODEGA(ndega).E_rel EPR_eval_MODEGA(ndega).E_rel EPR_out_MOGA(nga_MO).E_rel EPR_eval_MOGA(nga_MO).E_rel EPR_out_MODE(nde_MO).E_rel EPR_eval_MODE(nde_MO).E_rel; ...
    ];
fprintf(fileID,'%12s\r\n',' ');
fprintf(fileID,'%12s\r\n','Statistics of the optimal models');
fprintf(fileID,'%12s %12s %12s\n','MODEGA','MOGA','MODE');
fprintf(fileID,'%8s %8s %8s %8s %8s %8s\n','Training','Test','Training','Test','Training','Test');
fprintf(fileID,'%6.3f & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f\n',Tab_DE');

% compromise programming results
Tab_MO = [((1:DE.max)+1)' L_dega L_ga_MO];
fprintf(fileID,'%12s\r\n',' ');
fprintf(fileID,'%12s\r\n','Results of the multi-objective analyses:');
fprintf(fileID,'%d & %6.2f & %6.2f & %6.2f & %6.2f & & %6.2f & %6.2f & %6.2f & %6.2f\\\\\n',Tab_MO');
fclose(fileID);


for i = 1:DE.max
    % output file
    fileID = fopen('output_EPR.txt','a');
    fprintf(fileID,'%12s\r\n','---------- MO_DE_EPR Toolbox output file ----------');
    fprintf(fileID,'%12s\r\n','algorithm: DEGA');
    fprintf(fileID,'%12s %6.2f\r\n','EPR terms =',i+1);
    fprintf(fileID,'%12s\r\n','Best population');
    fprintf(fileID,'%6.2f %6.2f %6.2f \n',EPR_out_MODEGA(:,i).bestpop');
    fprintf(fileID,'%12s\r\n','Best parameters');
    fprintf(fileID,'%6.6f \n',EPR_out_MODEGA(:,i).best_par');
    fprintf(fileID,'%12s\r\n','Statistics');
    fprintf(fileID,'%12s %6.2f\r\n','R^2 =',EPR_out_MODEGA(:,i).R2);
    fprintf(fileID,'%12s %6.2f\r\n','RMSE =',EPR_out_MODEGA(:,i).RMSE);
    fprintf(fileID,'%12s %6.2f\r\n','MAE =',EPR_out_MODEGA(:,i).MAE);
    %fprintf(fileID,'%12s %6.2f\r\n','NSE =',EPR_out_DEGA(:,ndega).R2);
    fprintf(fileID,'%12s %6.2f\r\n','r =',EPR_out_MODEGA(:,i).r);
    %fprintf(fileID,'%12s %6.2f\r\n','SSE =',EPR_out_DEGA(:,ndega).sse_min);
    fprintf(fileID,'%12s %6.2f\r\n','E_rel =',EPR_out_MODEGA(:,i).E_rel);
    fclose(fileID);
end

    % output file of the best simulation
    fileID = fopen('output_EPR.txt','a');
    fprintf(fileID,'%12s\r\n','---------- MO_DE_EPR Toolbox output file ----------');
    fprintf(fileID,'%12s\r\n','algorithm: DEGA');
    fprintf(fileID,'%12s %6.2f\r\n','EPR terms =',ndega+1);
    fprintf(fileID,'%12s\r\n','Best population');
    fprintf(fileID,'%6.2f %6.2f %6.2f \n',EPR_out_MODEGA(:,ndega).bestpop);
    fprintf(fileID,'%12s\r\n','Best parameters');
    fprintf(fileID,'%6.6f \n',EPR_out_MODEGA(:,ndega).best_par');
    fprintf(fileID,'%12s\r\n','Statistics');
    fprintf(fileID,'%12s %6.2f\r\n','R^2 =',EPR_out_MODEGA(:,ndega).R2);
    fprintf(fileID,'%12s %6.2f\r\n','RMSE =',EPR_out_MODEGA(:,ndega).RMSE);
    fprintf(fileID,'%12s %6.2f\r\n','MAE =',EPR_out_MODEGA(:,ndega).MAE);
    %fprintf(fileID,'%12s %6.2f\r\n','NSE =',EPR_out_DEGA(:,ndega).R2);
    fprintf(fileID,'%12s %6.2f\r\n','r =',EPR_out_MODEGA(:,ndega).r);
    %fprintf(fileID,'%12s %6.2f\r\n','SSE =',EPR_out_DEGA(:,ndega).sse_min);
    fprintf(fileID,'%12s %6.2f\r\n','E_rel =',EPR_out_MODEGA(:,ndega).E_rel);
    fclose(fileID);
    

%end


close(hw);
elapsetime = toc; % stop timer
disp(['Elapse time (Multi-obj. EPR):       ' num2str(elapsetime)]); % display elapse time


%% Sub-plot
np = DE.max;
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
%--------------------- Plot A -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,2},...
'XGrid','off',...
'FontSize',fontsize,...
'Layer','top');
%% Post processing Compromise programming (R2)
% x-axis : number of EPR terms
ex = [1:np];
% x-axis : coefficient of determination (R2)
ey_MODE = A(1,1:np); % SODE
ey_MOGA = A(1,np+1:2*np); % SOGA
ey_DEGA = A(1,2*np+1:3*np); % MODE
%ey_GA_MO = A(1,3*np+1:4*np); % MOGA
%ey_DEGA = A(1,4*np+1:end-1); % DEGA
%ey_lit = A(1,end); % Published paper
% Now plot
fax = fontsize;
color_vec = [1,1,1];
% Create axes
% axes1 = axes(...
% 'Position',[0.169642857142857 0.184063494598109 0.769285714285714 0.752841267306653]);
% hold(axes1,'on');
% plot DEGA
P1 = plot(ex,ey_DEGA,'ks','LineStyle','-','Linewidth',1);
hold on
set(P1,'MarkerFaceColor',color_vec,...
'MarkerSize',8);
P2 = plot(ex,ey_MODE,'b^','LineStyle','-','Linewidth',1);
hold on
set(P2,'MarkerFaceColor',color_vec,...
'MarkerSize',8);
P3 = plot(ex,ey_MOGA,'ro','LineStyle','-','Linewidth',1);
hold on
set(P3,'MarkerFaceColor',color_vec,...
'MarkerSize',8);
% Plot published paper data
% P3 = plot(numel(Output_literatura.best_par)-1,ey_lit,'xk','Linewidth',2);
% hold on
% set(P3,'MarkerFaceColor','k',...
% 'MarkerSize',10);
P11 = plot(ex(ndega),ey_DEGA(ndega),'ks','LineStyle','-','Linewidth',1);
set(P11,'MarkerFaceColor','k',...
'MarkerSize',8);
P22 = plot(ex(nde_MO),ey_MODE(nde_MO),'b^','LineStyle','-','Linewidth',1);
set(P22,'MarkerFaceColor','b',...
'MarkerSize',8);
P33 = plot(ex(nga_MO),ey_MOGA(nga_MO),'ro','LineStyle','-','Linewidth',1);
set(P33,'MarkerFaceColor','r',...
'MarkerSize',8);
% plot grid
grid on
% axes adjustment
set(gca, ...
'Box' , 'on' , ...
'TickDir' , 'out' , ...
'TickLength' , [.02 .02] , ...
'XMinorTick' , 'on' , ...
'YMinorTick' , 'on' , ...
'LineWidth' , 1 );
% axes labels
xlabel('$m$','interpreter','latex','FontSize',fax);
ylabel('$R^2$ ','interpreter','latex','FontSize',fax);
set(gca, 'XTick', [1:np],'XTickLabel',{'1','2','3','4','5','6','7'});
set(gca,'fontsize',fax);
set(gca,'XLim',[1 DE.max]);
set(gca,'YLim',[0.7 1]);
set(gca,'TickLabelInterpreter','latex')
text(0.05,0.90,'(a)','units','normalized','FontSize',fax,'interpreter','latex');
%--------------------- Plot B -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,1},...
'XGrid','off',...
'FontSize',fontsize,...
'Layer','top');
% y-axis : coefficient of determination (R2)
ey_MODE = A(2,1:np); % SODE
ey_MOGA = A(2,np+1:2*np); % SOGA
ey_DEGA = A(2,2*np+1:3*np); % MODE
%ey_GA_MO = A(1,3*np+1:4*np); % MOGA
%ey_DEGA = A(1,4*np+1:end-1); % DEGA
%ey_lit = A(2,end); % Published paper
% Now plot
fax = fontsize;
color_vec = [1,1,1];
P1 = plot(ex,ey_DEGA,'ks','LineStyle','-','Linewidth',1);
hold on
set(P1,'MarkerFaceColor',color_vec,...
'MarkerSize',8);
P2 = plot(ex,ey_MODE,'b^','LineStyle','-','Linewidth',1);
hold on
set(P2,'MarkerFaceColor',color_vec,...
'MarkerSize',8);
P3 = plot(ex,ey_MOGA,'ro','LineStyle','-','Linewidth',1);
hold on
set(P3,'MarkerFaceColor',color_vec,...
'MarkerSize',8);
% Plot published paper data
% P3 = plot(numel(Output_literatura.best_par)-1,ey_lit,'xk','Linewidth',2);
% hold on
% set(P3,'MarkerFaceColor','k',...
% 'MarkerSize',10);
P11 = plot(ex(ndega),ey_DEGA(ndega),'ks','LineStyle','-','Linewidth',1);
set(P11,'MarkerFaceColor','k',...
'MarkerSize',8);
P22 = plot(ex(nde_MO),ey_MODE(nde_MO),'b^','LineStyle','-','Linewidth',1);
set(P22,'MarkerFaceColor','b',...
'MarkerSize',8);
P33 = plot(ex(nga_MO),ey_MOGA(nga_MO),'ro','LineStyle','-','Linewidth',1);
set(P33,'MarkerFaceColor','r',...
'MarkerSize',8);
% plot grid
grid on
% axes adjustment
set(gca, ...
'Box' , 'on' , ...
'TickDir' , 'out' , ...
'TickLength' , [.02 .02] , ...
'XMinorTick' , 'on' , ...
'YMinorTick' , 'on' , ...
'LineWidth' , 1 );
% axes labels
xlabel('$m$','interpreter','latex','FontSize',fax);
ylabel('RMSE ','interpreter','latex','FontSize',fax);
set(gca, 'XTick', [1:np],'XTickLabel',{'1','2','3','4','5','6','7'});
set(gca,'fontsize',fax);
set(gca,'TickLabelInterpreter','latex')
set(gca,'XLim',[1 DE.max]);
%set(gca,'YLim',[0.7 1]);
text(0.05,0.90,'(b)','units','normalized','FontSize',fax,'interpreter','latex');
%--------------------- Plot C -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,2},...
'XGrid','off',...
'FontSize',fontsize,...
'Layer','top');
% y-axis : coefficient of determination (R2)
ey_MODE = A(3,1:np); % SODE
ey_MOGA = A(3,np+1:2*np); % SOGA
ey_DEGA = A(3,2*np+1:3*np); % MODE
%ey_GA_MO = A(1,3*np+1:4*np); % MOGA
%ey_DEGA = A(1,4*np+1:end-1); % DEGA
%ey_lit = A(3,end); % Published paper
% Now plot
fax = fontsize;
color_vec = [1,1,1];
% plot DEGA
P1 = plot(ex,ey_DEGA,'ks','LineStyle','-','Linewidth',1);
hold on
set(P1,'MarkerFaceColor',color_vec,...
'MarkerSize',8);
P2 = plot(ex,ey_MODE,'b^','LineStyle','-','Linewidth',1);
hold on
set(P2,'MarkerFaceColor',color_vec,...
'MarkerSize',8);
P3 = plot(ex,ey_MOGA,'ro','LineStyle','-','Linewidth',1);
hold on
set(P3,'MarkerFaceColor',color_vec,...
'MarkerSize',8);
% Plot published paper data
% P3 = plot(numel(Output_literatura.best_par)-1,ey_lit,'xk','Linewidth',2);
% hold on
% set(P3,'MarkerFaceColor','k',...
% 'MarkerSize',10);
P11 = plot(ex(ndega),ey_DEGA(ndega),'ks','LineStyle','-','Linewidth',1);
set(P11,'MarkerFaceColor','k',...
'MarkerSize',8);
P22 = plot(ex(nde_MO),ey_MODE(nde_MO),'b^','LineStyle','-','Linewidth',1);
set(P22,'MarkerFaceColor','b',...
'MarkerSize',8);
P33 = plot(ex(nga_MO),ey_MOGA(nga_MO),'ro','LineStyle','-','Linewidth',1);
set(P33,'MarkerFaceColor','r',...
'MarkerSize',8);
% plot grid
grid on
% axes adjustment
set(gca, ...
'Box' , 'on' , ...
'TickDir' , 'out' , ...
'TickLength' , [.02 .02] , ...
'XMinorTick' , 'on' , ...
'YMinorTick' , 'on' , ...
'LineWidth' , 1 );
% axes labels
xlabel('$m$','interpreter','latex','FontSize',fax);
ylabel('$r$ ','interpreter','latex','FontSize',fax);
set(gca, 'XTick', [1:np],'XTickLabel',{'1','2','3','4','5','6','7'});
set(gca,'fontsize',fax);
set(gca,'TickLabelInterpreter','latex')
text(0.05,0.90,'(c)','units','normalized','FontSize',fax,'interpreter','latex');
set(gca,'XLim',[1 DE.max]);
%--------------------- Plot D -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,1},...
'XGrid','off',...
'FontSize',fontsize,...
'Layer','top');
%% Post processing Compromise programming (E_rel)
% y-axis : coefficient of determination (R2)
ey_MODE = A(4,1:np); % SODE
ey_MOGA = A(4,np+1:2*np); % SOGA
ey_DEGA = A(4,2*np+1:3*np); % MODE
%ey_GA_MO = A(1,3*np+1:4*np); % MOGA
%ey_DEGA = A(1,4*np+1:end-1); % DEGA
%ey_lit= A(4,end); % Published paper
% Now plot
fax = fontsize;
color_vec = [1,1,1];
% plot DEGA
P1 = plot(ex,ey_DEGA,'ks','LineStyle','-','Linewidth',1);
hold on
set(P1,'MarkerFaceColor',color_vec,...
'MarkerSize',8);
P2 = plot(ex,ey_MODE,'b^','LineStyle','-','Linewidth',1);
hold on
set(P2,'MarkerFaceColor',color_vec,...
'MarkerSize',8);
P3 = plot(ex,ey_MOGA,'ro','LineStyle','-','Linewidth',1);
hold on
set(P3,'MarkerFaceColor',color_vec,...
'MarkerSize',8);
% Plot published paper data
% P3 = plot(numel(Output_literatura.best_par)-1,ey_lit,'xk','Linewidth',2);
% hold on
% set(P3,'MarkerFaceColor','k',...
% 'MarkerSize',10);
P11 = plot(ex(ndega),ey_DEGA(ndega),'ks','LineStyle','-','Linewidth',1);
set(P11,'MarkerFaceColor','k',...
'MarkerSize',8);
P22 = plot(ex(nde_MO),ey_MODE(nde_MO),'b^','LineStyle','-','Linewidth',1);
set(P22,'MarkerFaceColor','b',...
'MarkerSize',8);
P33 = plot(ex(nga_MO),ey_MOGA(nga_MO),'ro','LineStyle','-','Linewidth',1);
set(P33,'MarkerFaceColor','r',...
'MarkerSize',8);
% plot grid
grid on
% axes adjustment
set(gca, ...
'Box' , 'on' , ...
'TickDir' , 'out' , ...
'TickLength' , [.02 .02] , ...
'XMinorTick' , 'on' , ...
'YMinorTick' , 'on' , ...
'LineWidth' , 1 );
% axes labels
xlabel('$m$','interpreter','latex','FontSize',fax);
ylabel('$E_{\rm rel}$ ','interpreter','latex','FontSize',fax);
set(gca, 'XTick', [1:np],'XTickLabel',{'1','2','3','4','5','6','7'});
set(gca,'XLim',[1 DE.max]);
set(gca,'fontsize',fax);
set(gca,'TickLabelInterpreter','latex')
text(0.05,0.90,'(d)','units','normalized','FontSize',fax,'interpreter','latex');

%% Sub-plot

% parameters for figure and panel size
plotheight=12;
plotwidth=8;
subplotsx=1;
subplotsy=3;   
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
ax=axes('position',sub_pos{1,3},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');

%% Post processing EPR (SSE)
sv = 1; % space vector

fax = fontsize;
color_vec = [1,1,1];

P3 = plot((1:sv:DE.last)',EPR_out_MODE(:,nde_MO).SSE(1:sv:end),'color',[.5 .5 .5],'LineStyle','--','Linewidth',2);  
hold on

P4 = plot((1:sv:DE.last)',EPR_out_MOGA(:,nga_MO).SSE(1:sv:end),'r-','LineStyle','-.','Linewidth',2);  

P5 = plot((1:sv:DE.last)',EPR_out_MODEGA(:,ndega).SSE(1:sv:end),'k-','LineStyle','-','Linewidth',2);  

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
   
%set(gca,'YLim',[5 15]);
set(gca,'xticklabel',{[]});
leg1 = legend([ P5, P4, P3],['MODEGA'],['MOGA'],['MODE'],...
'location', 'best','Orientation','horizontal');
set(leg1,'Interpreter','latex');
text(0.9,0.85,'(a)','units','normalized','FontSize',fax,'interpreter','latex');   

%set(gcf, 'Units', 'pixels', 'Position', [10, 100, 1000, 400]);
%--------------------- Plot C -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');

%% Post processing EPR (Number of variables)

MOGA = EPR_out_MOGA(:,nga_MO).OF(:,2);
MODE = EPR_out_MODE(:,nde_MO).OF(:,2);
MODEGA  = EPR_out_MODEGA(:,ndega).OF(:,2);


sv = 1; % space vector

%fax = 18;
color_vec = [1,1,1];

P3 = plot((1:sv:DE.last)',MODE,'color',[.5 .5 .5],'LineStyle','--','Linewidth',2);  
hold on

P4 = plot((1:sv:DE.last)',MOGA,'r-','LineStyle','-.','Linewidth',2);  

P5 = plot((1:sv:DE.last)',MODEGA,'k-','LineStyle','-','Linewidth',2);  

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
%xlabel('Number of generations ','interpreter','latex','FontSize',fax);
ylabel({'Number of'; 'variables'},'interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);
set(leg1,'Interpreter','latex');    
set(gca,'YLim',[2 7]);
set(gca,'xticklabel',{[]});
text(0.9,0.85,'(b)','units','normalized','FontSize',fax,'interpreter','latex');   
%--------------------- Plot C -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');

%% Post processing EPR (self adaptive offspring creation)

GA_plot = EPR_out_MODEGA(:,ndega).N_alg(:,2);
DE_plot = EPR_out_MODEGA(:,ndega).N_alg(:,1);

color_vec = [1,1,1];
% Create axes
sv = 1; % space vector

%fax = 18;
color_vec = [1,1,1];
% Create axes

P3 = plot((1:sv:DE.last)',DE_plot(1:sv:end),'b','LineStyle','-','Linewidth',1);  
hold on
% set(P3,'MarkerFaceColor',color_vec,...
%     'MarkerSize',8)
hold on
P4 = plot((1:sv:DE.last)',GA_plot(1:sv:end),'r','LineStyle','--','Linewidth',1);  
hold on

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
xlabel('Number of generations','interpreter','latex','FontSize',fax);
ylabel('$N_o$','interpreter','latex','FontSize',fax);
grid on
set(gca,'YLim',[0 60]);
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);
leg1 = legend([P3, P4],['DE'],['GA'],...
'location', 'north','Orientation','horizontal');
set(leg1,'Interpreter','latex'); 
text(0.9,0.85,'(c)','units','normalized','FontSize',fax,'interpreter','latex');   


%% Sub-plot R2

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

P3 = plot(Y_cal,EPR_out_MODEGA(:,ndega).besty,'ks');
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

P5 = plot(Y_eval,EPR_eval_MODEGA(:,ndega).Y_eval,'ks');
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


% print final equation
fprintf('----------------------------------------------- \n')
fprintf('Equation (MODEGA) \n')
fprintf('\n')
print_EPR_equation (X,EPR_out_MODEGA(ndega).best_par',EPR_out_MODEGA(ndega).bestpop);        
fprintf('----------------------------------------------- \n')


% If Sensitivity-Driven is not selected then stop code
if  strcmp(MODEGASD, 'No')
    return
end

%% RUN CROSS-VALIDATION 


% Update number of EPR terms with optimum value (excluding bias) for MODEGA
DE.m = ndega; GA.m = ndega;         
DE.N = ndega*DE.pop;      % Population size

% Loop to create different random samplings
for i = 1:n_split
    
    if i == 1 % if first batch maintain article data
        EPR_out_CrossVal(i).Y_cal = Y_cal; % Store data
        EPR_out_CrossVal(i).Y_eval = Y_eval;
        EPR_out_CrossVal(i).X_cal  = X_cal;
        EPR_out_CrossVal(i).X_eval = X_eval;
        
    else    % Resample
        [m,n] = size(Y) ; % get lines (m) and column (n) sizes
        idx = randperm(m)  ; % shuffle indices

        % Get first P % and assign to calibration
        EPR_out_CrossVal(i).Y_cal(:,1) = Y(idx(1:round(P*m)),:) ; 
        % Get last P % and assign to evaluation
        EPR_out_CrossVal(i).Y_eval(:,1) = Y(idx(round(P*m)+1:end),:) ;

        % Loop to store independent variables into 3D matrix
        for j=1:size(X,2)
            EPR_out_CrossVal(i).X_cal(:,j) = X(idx(1:round(P*m)),j) ; 
            EPR_out_CrossVal(i).X_eval(:,j) = X(idx(round(P*m)+1:end),j) ;
        end

        % Get current input matrix
        Y_cal  = EPR_out_CrossVal(i).Y_cal;
        Y_eval = EPR_out_CrossVal(i).Y_eval;

        X_cal  = EPR_out_CrossVal(i).X_cal;
        X_eval = EPR_out_CrossVal(i).X_eval;
    end
    
    %% RUN MONTE CARLO
    
    %-------------------------------------------------------------------------%
    %-------------------------------------------------------------------------%
    % ------------------------ Monte Carlo Simulations -----------------------%
    %-------------------------------------------------------------------------%
    %-------------------------------------------------------------------------%

    % Here we run multiple EPR's with optimum m parameter based on compromise
    % programming Choice

    hw = waitbar(0,'Running MonteCarlo EPR-MODEGA...');

    tic % start timer

    % Run MODEGA Monte Carlo
    for j=1:n_MonteCarlo
            % update wait bar
           waitbar(j/n_MonteCarlo);

           %Run EPR-MODEGA
           %[Output_MonteCarlo_MODEGA] = MODEGA_algorithm (DE, GA, Y_cal, X_cal);
           [Output_MonteCarlo_MODEGA] = MODEGA_algorithm_Cal_Eval(DE, GA, Y_cal, X_cal, Y_eval, X_eval);

           % save data in a structure variable
            EPR_out_MonteCarlo_MODEGA(:,j).SSE = Output_MonteCarlo_MODEGA.SEE;           % best SSE for each generation 
            EPR_out_MonteCarlo_MODEGA(:,j).besty = Output_MonteCarlo_MODEGA.besty;       % best simulation
            EPR_out_MonteCarlo_MODEGA(:,j).bestpop = Output_MonteCarlo_MODEGA.bestpop;   % best population
            EPR_out_MonteCarlo_MODEGA(:,j).best_par = Output_MonteCarlo_MODEGA.best_par; % best EPR parameters
            EPR_out_MonteCarlo_MODEGA(:,j).R2 = Output_MonteCarlo_MODEGA.R2;             % R2
            EPR_out_MonteCarlo_MODEGA(:,j).RMSE = Output_MonteCarlo_MODEGA.RMSE;         % RMSE
            EPR_out_MonteCarlo_MODEGA(:,j).MAE = Output_MonteCarlo_MODEGA.MAE;           % MAE
            EPR_out_MonteCarlo_MODEGA(:,j).E_rel = Output_MonteCarlo_MODEGA.E_rel;       % Relative efficiency
            EPR_out_MonteCarlo_MODEGA(:,j).r = Output_MonteCarlo_MODEGA.r;               % Pearson's correlation 
            EPR_out_MonteCarlo_MODEGA(:,j).PBIAS = Output_MonteCarlo_MODEGA.PBIAS;       % Percentage Bias 
            EPR_out_MonteCarlo_MODEGA(:,j).OF = Output_MonteCarlo_MODEGA.OF;             % best OF-values for each generation (SSE and num. variables in the model)
            EPR_out_MonteCarlo_MODEGA(:,j).N_alg = Output_MonteCarlo_MODEGA.N_alg;       % number of offspring

            % save calibration statistics for each generation
            EPR_out_MonteCarlo_MODEGA(:,j).stats_cal(:,1) = Output_MonteCarlo_MODEGA.stats(:,1);       %R2
            EPR_out_MonteCarlo_MODEGA(:,j).stats_cal(:,2) = Output_MonteCarlo_MODEGA.stats(:,2);       %RMSE
            EPR_out_MonteCarlo_MODEGA(:,j).stats_cal(:,3) = Output_MonteCarlo_MODEGA.stats(:,3);   %MAE
            EPR_out_MonteCarlo_MODEGA(:,j).stats_cal(:,4) = Output_MonteCarlo_MODEGA.stats(:,4); %R_rel
            EPR_out_MonteCarlo_MODEGA(:,j).stats_cal(:,5) = Output_MonteCarlo_MODEGA.stats(:,5);     %r
            EPR_out_MonteCarlo_MODEGA(:,j).stats_cal(:,6) = Output_MonteCarlo_MODEGA.stats(:,6);     %PBIAS
            
            % save evaluation statistics for each generation
            EPR_out_MonteCarlo_MODEGA(:,j).stats_eval(:,1) = Output_MonteCarlo_MODEGA.stats_eval(:,1);       %R2
            EPR_out_MonteCarlo_MODEGA(:,j).stats_eval(:,2) = Output_MonteCarlo_MODEGA.stats_eval(:,2);       %RMSE
            EPR_out_MonteCarlo_MODEGA(:,j).stats_eval(:,3) = Output_MonteCarlo_MODEGA.stats_eval(:,3);   %MAE
            EPR_out_MonteCarlo_MODEGA(:,j).stats_eval(:,4) = Output_MonteCarlo_MODEGA.stats_eval(:,4); %R_rel
            EPR_out_MonteCarlo_MODEGA(:,j).stats_eval(:,5) = Output_MonteCarlo_MODEGA.stats_eval(:,5);     %r
            EPR_out_MonteCarlo_MODEGA(:,j).stats_eval(:,6) = Output_MonteCarlo_MODEGA.stats_eval(:,6);     %PBIAS

            % Run EPR equation with evaluation data set
            [Y_MonteCarlo_MODEGA_eval] = evaluation_DE_regression (DE, Output_MonteCarlo_MODEGA, X_eval);
            EPR_eval_MonteCarlo_MODEGA(:,j).Y_eval = Y_MonteCarlo_MODEGA_eval;
            % Compute evaluation statistics (best model)
            [Out_stat] = statistic_model (Y_eval, Y_MonteCarlo_MODEGA_eval);
            EPR_eval_MonteCarlo_MODEGA(:,j).R2 = Out_stat.R2;
            EPR_eval_MonteCarlo_MODEGA(:,j).RMSE = Out_stat.RMSE;
            EPR_eval_MonteCarlo_MODEGA(:,j).MAE = Out_stat.MAE;
            EPR_eval_MonteCarlo_MODEGA(:,j).E_rel = Out_stat.E_rel;
            EPR_eval_MonteCarlo_MODEGA(:,j).r = Out_stat.r;  


            % Run sensitivity and store output
            [X_sensitivity, Y_sensitivity] = sensitivity(DE, Output_MonteCarlo_MODEGA,X_cal);
            EPR_out_MonteCarlo_MODEGA(:,j).Y_sensitivity = Y_sensitivity;
         
            % ------------------ CROSS VALIDATION DATA -------------------%
            % Save for Cross-validation - First line gen 1 Last line gen 300
            EPR_out_CrossVal(i).SSE(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).SSE;
            % Save DE
            EPR_out_CrossVal(i).N_alg_DE(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).N_alg(:,1);            
            % Save GA
            EPR_out_CrossVal(i).N_alg_GA(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).N_alg(:,2); 
            
             % Save Monte Carlo Statistics
            % save calibration statistics for each generation
            EPR_out_CrossVal(i).Gen_Cal_R2(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).stats_cal(:,1);      %R2
            EPR_out_CrossVal(i).Gen_Cal_RMSE(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).stats_cal(:,2);    %RMSE
            EPR_out_CrossVal(i).Gen_Cal_MAE(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).stats_cal(:,3);     %MAE
            EPR_out_CrossVal(i).Gen_Cal_E_rel(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).stats_cal(:,4);   %E_rel
            EPR_out_CrossVal(i).Gen_Cal_r(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).stats_cal(:,5);       %r
            EPR_out_CrossVal(i).Gen_Cal_PBIAS(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).stats_cal(:,6);   %PBIAS

            % save evaluation statistics for each generation
            EPR_out_CrossVal(i).Gen_eval_R2(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).stats_eval(:,1);     %R2
            EPR_out_CrossVal(i).Gen_eval_RMSE(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).stats_eval(:,2);     %RMSE
            EPR_out_CrossVal(i).Gen_eval_MAE(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).stats_eval(:,3);     %MAE
            EPR_out_CrossVal(i).Gen_eval_E_rel(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).stats_eval(:,4);     %E_rel
            EPR_out_CrossVal(i).Gen_eval_r(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).stats_eval(:,5);     %r
            EPR_out_CrossVal(i).Gen_eval_PBIAS(:,j) = EPR_out_MonteCarlo_MODEGA(:,j).stats_eval(:,6);     %PBIAS

            % save sensitivity
            EPR_out_CrossVal(i).Y_sens_MonteCarlo(:,:,j) = EPR_out_MonteCarlo_MODEGA(j).Y_sensitivity;
    % end of Monte Carlo        
    end
    
    % Get best model info from each cross-validation
    % Save best model of Monte Carlo in struct variable
    [BEST] = model_selection(n_MonteCarlo,X_sensitivity,EPR_out_MonteCarlo_MODEGA, X_cal, Y_cal, Names);
    
    % Save best equation and its statistics
    EPR_out_CrossVal(i).besty = EPR_out_MonteCarlo_MODEGA(:,BEST).besty  ;      
    EPR_out_CrossVal(i).bestpop = EPR_out_MonteCarlo_MODEGA(:,BEST).bestpop  ;
    EPR_out_CrossVal(i).best_par = EPR_out_MonteCarlo_MODEGA(:,BEST).best_par  ;
    EPR_out_CrossVal(i).R2 = EPR_out_MonteCarlo_MODEGA(:,BEST).R2  ;    
    EPR_out_CrossVal(i).RMSE = EPR_out_MonteCarlo_MODEGA(:,BEST).RMSE  ;  
    EPR_out_CrossVal(i).r = EPR_out_MonteCarlo_MODEGA(:,BEST).r  ;
    EPR_out_CrossVal(i).E_rel = EPR_out_MonteCarlo_MODEGA(:,BEST).E_rel  ;
    EPR_out_CrossVal(i).PBIAS = EPR_out_MonteCarlo_MODEGA(:,BEST).PBIAS  ;     
    % Save best sensitivity
    EPR_out_CrossVal(i).Y_sens_SelectModel = EPR_out_MonteCarlo_MODEGA(:,BEST).Y_sensitivity;
    EPR_out_CrossVal(i).BEST = BEST;
    
    close(hw);
    elapsetime = toc; % stop timer
    disp(['Elapse time (MonteCarlo):       ' num2str(elapsetime)]); % display elapse time
    
    % If Cross validation is not selected the plot results and stop code
    if n_split == 1 || strcmp(Data_divided, 'Yes') && i == 1
        % plot uncertainty graphs
        %[p]= plot_MonteCarlo_uncertainty(ndega, DE.last, n_MonteCarlo, EPR_out_MODEGA,EPR_out_MonteCarlo_MODEGA,EPR_eval_MonteCarlo_MODEGA);
        [p]= plot_MonteCarlo_uncertainty(DE.last, n_MonteCarlo, EPR_out_MonteCarlo_MODEGA,EPR_eval_MonteCarlo_MODEGA,BEST, Y_cal, Y_eval, Names);
        % Plot sensitivity graphs
        [mean_sensitivity]=plot_sensitivity_ranges(n_MonteCarlo, BEST ,X_sensitivity,EPR_out_MonteCarlo_MODEGA, Names);
        
        % output file of the best simulation
        fileID = fopen('MonteCarlo_EPR.txt','a');
        fprintf(fileID,'%12s\r\n','---------- MonteCarlo MODEGA Toolbox output file ----------');
        fprintf(fileID,'%12s\r\n','--------------------------------------');
        fprintf(fileID,'%12s\r\n','algorithm: DEGA');
        fprintf(fileID,'%12s %6.2f\r\n','EPR terms (+ bias) =',ndega+1);
        fprintf(fileID,'%12s %6.2f\r\n','Best model No. =',BEST);
        fprintf(fileID,'%12s\r\n','Best population');
        for iii=1:ndega
            fprintf(fileID, [repmat('%6.2f',1,DE.k), '\n'],EPR_out_MonteCarlo_MODEGA(BEST).bestpop(iii,:));
        end
        fprintf(fileID,'%12s\r\n','Best parameters');
        fprintf(fileID,'%6.6f \n',EPR_out_MonteCarlo_MODEGA(BEST).best_par');
        fprintf(fileID,'%12s\r\n','--------------------------------------');
        fprintf(fileID,'%12s\r\n','Calibration statistics');
        fprintf(fileID,'%12s %6.3f\r\n','R^2 =',EPR_out_MonteCarlo_MODEGA(BEST).R2);
        fprintf(fileID,'%12s %6.3f\r\n','RMSE =',EPR_out_MonteCarlo_MODEGA(BEST).RMSE);
        fprintf(fileID,'%12s %6.3f\r\n','MAE =',EPR_out_MonteCarlo_MODEGA(BEST).MAE);
        fprintf(fileID,'%12s %6.3f\r\n','r =',EPR_out_MonteCarlo_MODEGA(BEST).r);
        fprintf(fileID,'%12s %6.3f\r\n','E_rel =',EPR_out_MonteCarlo_MODEGA(BEST).E_rel);
        fprintf(fileID,'%12s\r\n','--------------------------------------');
        fprintf(fileID,'%12s\r\n','Evaluation statistics');
        fprintf(fileID,'%12s %6.3f\r\n','R^2 =',EPR_eval_MonteCarlo_MODEGA(BEST).R2);
        fprintf(fileID,'%12s %6.3f\r\n','RMSE =',EPR_eval_MonteCarlo_MODEGA(BEST).RMSE);
        fprintf(fileID,'%12s %6.3f\r\n','MAE =',EPR_eval_MonteCarlo_MODEGA(BEST).MAE);
        fprintf(fileID,'%12s %6.3f\r\n','r =',EPR_eval_MonteCarlo_MODEGA(BEST).r);
        fprintf(fileID,'%12s %6.3f\r\n','E_rel =',EPR_eval_MonteCarlo_MODEGA(BEST).E_rel);
        fprintf(fileID,'%12s\r\n','--------------------------------------');
        fclose(fileID);
        
        % print final equation
        fprintf('----------------------------------------------- \n')
        fprintf('Equation (MODEGA-SD) \n')
        fprintf('\n')
        print_EPR_equation (X,EPR_out_MonteCarlo_MODEGA(BEST).best_par',EPR_out_MonteCarlo_MODEGA(BEST).bestpop);        
        fprintf('----------------------------------------------- \n')
        
        % Stop code
        if strcmp(Data_divided, 'No')
            return
        end
    end
        
    % Else move into different training and testing combination   
end

%% PLOT

plotheight=12;%12
plotwidth=10;%8
subplotsx=1;
subplotsy=4;   
leftedge=1.5;
rightedge=0.5;   
topedge=0.5;
bottomedge=1.5;
spacex=1.2;
spacey=1.4;
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

% Calcula valor médio de SSE para primeira combinacao de dados
Media = mean(EPR_out_CrossVal(1).SSE,2);
% Obter intervalos de confiança (95%)
[SSE.low, SSE.up] = confidence_intervals (EPR_out_CrossVal(1).SSE');

% Plotar Incerteza (errorbar é utilizado apenas para legenda)
n_gen_vec=linspace(1,DE.last,DE.last);
P3 = Fill_Ranges(n_gen_vec,SSE.low,SSE.up, [0. 0. 0.]); hold on
P4 = plot((1:1:DE.last)',Media,'r-','LineStyle','-','Linewidth',2,'Color',[0. 0. 0.]); hold on;
% Error bar for legend
error1 = errorbar(-20,1,1,'o-','Color',[0. 0. 0.],...
                    'MarkerFaceColor','w','LineWidth',1.5,'LineStyle','none'); hold on;

% Plotagem manual da media e intervalo de confiança para 50 gerações
plot([50-3, 50+3],[SSE.up(1,50),SSE.up(1,50)],'k-','Linewidth',1.5,'Color',[0. 0. 0.]);hold on;
plot([50-3, 50+3],[SSE.low(1,50),SSE.low(1,50)],'k-','Linewidth',1.5,'Color',[0. 0. 0.]);hold on;
plot([50,50],[Media(50,1), SSE.low(1,50)],'k-','Linewidth',1.5,'Color',[0. 0. 0.]);hold on;
plot([50,50],[Media(50,1), SSE.up(1,50)],'k-','Linewidth',1.5,'Color',[0. 0. 0.]);hold on;
plot(50,Media(50,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',[0. 0. 0.]);hold on;

% Plotagem manual da media e intervalo de confiança para 200 gerações
plot([200-3, 200+3],[SSE.up(1,200),SSE.up(1,200)],'k-','Linewidth',1.5,'Color',[0. 0. 0.]);hold on;
plot([200-3, 200+3],[SSE.low(1,200),SSE.low(1,200)],'k-','Linewidth',1.5,'Color',[0. 0. 0.]);hold on;
plot([200,200],[Media(200,1), SSE.low(1,200)],'k-','Linewidth',1.5,'Color',[0. 0. 0.]);hold on;
plot([200,200],[Media(200,1), SSE.up(1,200)],'k-','Linewidth',1.5,'Color',[0. 0 0]);hold on;
plot(200,Media(200,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',[0 0 0]);hold on;

% Atribui transparencia para intervalo de confiança
alpha(P3,.3);

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );

xlabel('Number of generations ','interpreter','latex','FontSize',fax);
ylabel('SSR ($\%^2$)','interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);
   
set(gca,'YLim',[5 15]);
set(gca,'XLim',[0 300]);
%set(gca,'xticklabel',{[]});
leg1 = legend([P4, P3 , error1],['Mean'],['95$\%$ uncertainty'],['95$\%$ uncertainty'],...
'location', 'north','Orientation','horizontal');
set(leg1,'Interpreter','latex');
text(0.95,0.85,'(a)','units','normalized','FontSize',fax,'interpreter','latex');   


%--------------------- Plot B -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,3},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');

for i = 1:n_split
    Media = mean(EPR_out_CrossVal(i).SSE,2);
    [SSE.low, SSE.up] = confidence_intervals (EPR_out_CrossVal(i).SSE');
    plot([i-0.12, i+0.12],[SSE.up(1,50),SSE.up(1,50)],'k-','Linewidth',1.5,'Color',[0 0 0]);hold on;
    plot([i-0.12, i+0.12],[SSE.low(1,50),SSE.low(1,50)],'k-','Linewidth',1.5,'Color',[0 0 0]);hold on;
    plot([i,i],[Media(50,1), SSE.low(1,50)],'k-','Linewidth',1.5,'Color',[0 0 0]);hold on;
    plot([i,i],[Media(50,1), SSE.up(1,50)],'k-','Linewidth',1.5,'Color',[0 0 0]);hold on;
    plot(i,Media(50,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',[0 0 0]);hold on;
end

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
xlabel('Cross-validation batch','interpreter','latex','FontSize',fax);
ylabel('SSR ($\%^2$)','interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);

set(gca,'XLim',[0.5 n_split+0.5]);
text(0.95,0.85,'(b)','units','normalized','FontSize',fax,'interpreter','latex');  
%text(0.1,0.85,'Generation 50','units','normalized','FontSize',fax,'interpreter','latex');  


%--------------------- Plot C -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');

for i = 1:n_split
    Media = mean(EPR_out_CrossVal(i).SSE,2);
    [SSE.low, SSE.up] = confidence_intervals (EPR_out_CrossVal(i).SSE');
    plot([i-0.12, i+0.12],[SSE.up(1,200),SSE.up(1,200)],'k-','Linewidth',1.5,'Color',[0 0 0]);hold on;
    plot([i-0.12, i+0.12],[SSE.low(1,200),SSE.low(1,200)],'k-','Linewidth',1.5,'Color',[0 0 0]);hold on;
    plot([i,i],[Media(200,1), SSE.low(1,200)],'k-','Linewidth',1.5,'Color',[0 0 0]);hold on;
    plot([i,i],[Media(200,1), SSE.up(1,200)],'k-','Linewidth',1.5,'Color',[0 0 0]);hold on;
    plot(i,Media(200,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',[0 0 0]);hold on;
end

set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
xlabel('Cross-validation batch','interpreter','latex','FontSize',fax);
ylabel('SSR ($\%^2$)','interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);

set(gca,'XLim',[0.5 n_split+0.5]);
text(0.95,0.85,'(c)','units','normalized','FontSize',fax,'interpreter','latex');  
%text(0.1,0.85,'Generation 200','units','normalized','FontSize',fax,'interpreter','latex');  



%--------------------- Plot D -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');

for i = 1:n_split
   Media_DE = mean(EPR_out_CrossVal(i).N_alg_DE,2);
   P4 = plot((1:1:DE.last)',Media_DE,'b-','LineStyle','-','Linewidth',1.5,'Color',[0.7 0.2 1]);hold on;
   Media_GA = mean(EPR_out_CrossVal(i).N_alg_GA,2);
   P5 = plot((1:1:DE.last)',Media_GA,'r-','LineStyle','-','Linewidth',1.5,'Color',[1 0.4 .2]);hold on;
end


      
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
xlabel('Generations','interpreter','latex','FontSize',fax);
ylabel('N$_{\rm o}$','interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);

set(gca,'XLim',[0 DE.last]);
set(gca,'YLim',[0 50]);

leg1 = legend([P4, P5 ],['DE mean'],['GA mean'],...
'location', 'north','Orientation','horizontal');
set(leg1,'Interpreter','latex');

text(0.95,0.85,'(d)','units','normalized','FontSize',fax,'interpreter','latex'); 



%-------------------------- Train - Test ---------------------------------%



plotheight=5;%12
plotwidth=17;%8
subplotsx=2;
subplotsy=2;   
leftedge=1.5;
rightedge=0.5;   
topedge=0.5;
bottomedge=1.5;
spacex=1.7;
spacey=0.6;
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
ax=axes('position',sub_pos{1,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');

Color_cal  = [0.2 0.7 1];
Color_eval = [1 0.7 0.2];


% Calibration
for i = 1:n_split
    % Calcula valor médio de SSE para primeira combinacao de dados
    Media = mean(EPR_out_CrossVal(i).Gen_Cal_R2,2);
    % Obter intervalos de confiança (95%)
    [SSE.low, SSE.up] = confidence_intervals (EPR_out_CrossVal(i).Gen_Cal_R2');
    
    j=i-0.15;
    plot([j-0.12, j+0.12],[SSE.up(1,DE.last),SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot([j-0.12, j+0.12],[SSE.low(1,DE.last),SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot([j,j],[Media(DE.last,1), SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot([j,j],[Media(DE.last,1), SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot(j,Media(DE.last,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_cal);hold on;
end


% Eval
for i = 1:n_split
    % Calcula valor médio de SSE para primeira combinacao de dados
    Media = mean(EPR_out_CrossVal(i).Gen_eval_R2,2);
    % Obter intervalos de confiança (95%)
    [SSE.low, SSE.up] = confidence_intervals (EPR_out_CrossVal(i).Gen_eval_R2');
    
    j=i+0.15;
    plot([j-0.12, j+0.12],[SSE.up(1,DE.last),SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot([j-0.12, j+0.12],[SSE.low(1,DE.last),SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot([j,j],[Media(DE.last,1), SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot([j,j],[Media(DE.last,1), SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot(j,Media(DE.last,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_eval);hold on;
end

% Error bar for legend
error1 = errorbar(-20,1,1,'o-','Color',Color_cal,...
                    'MarkerFaceColor','w','LineWidth',1.5,'LineStyle','none'); hold on;
% Error bar for legend
error2 = errorbar(-20,1,1,'o-','Color',Color_eval,...
                    'MarkerFaceColor','w','LineWidth',1.5,'LineStyle','none'); hold on;


set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
xlabel('Cross-validation batch','interpreter','latex','FontSize',fax);
ylabel('R$^2$','interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);

leg1 = legend([error1, error2],['Train'],['Test'],...
'location', 'southwest','Orientation','vertical');
set(leg1,'Interpreter','latex');


set(gca,'XLim',[1-0.5 i+0.5]);
xticks([1 2 3 4 5 6 7 8 9 10])
text(0.905,0.925,'(a)','units','normalized','FontSize',fax,'interpreter','latex'); 


%--------------------- Plot B -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,2},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');

% Calibration
for i = 1:n_split
    % Calcula valor médio de SSE para primeira combinacao de dados
    Media = mean(EPR_out_CrossVal(i).Gen_Cal_RMSE,2);
    % Obter intervalos de confiança (95%)
    [SSE.low, SSE.up] = confidence_intervals (EPR_out_CrossVal(i).Gen_Cal_RMSE');
    
    j=i-0.15;
    plot([j-0.12, j+0.12],[SSE.up(1,DE.last),SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot([j-0.12, j+0.12],[SSE.low(1,DE.last),SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot([j,j],[Media(DE.last,1), SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot([j,j],[Media(DE.last,1), SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot(j,Media(DE.last,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_cal);hold on;
end

% Eval
for i = 1:n_split
    % Calcula valor médio de SSE para primeira combinacao de dados
    Media = mean(EPR_out_CrossVal(i).Gen_eval_RMSE,2);
    % Obter intervalos de confiança (95%)
    [SSE.low, SSE.up] = confidence_intervals (EPR_out_CrossVal(i).Gen_eval_RMSE');
    
    j=i+0.15;
    plot([j-0.12, j+0.12],[SSE.up(1,DE.last),SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot([j-0.12, j+0.12],[SSE.low(1,DE.last),SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot([j,j],[Media(DE.last,1), SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot([j,j],[Media(DE.last,1), SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot(j,Media(DE.last,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_eval);hold on;
end


set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
xlabel('Cross-validation batch','interpreter','latex','FontSize',fax);
ylabel('RMSE','interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);


set(gca,'XLim',[1-0.5 i+0.5]);

xticks([1 2 3 4 5 6 7 8 9 10])
text(0.905,0.925,'(b)','units','normalized','FontSize',fax,'interpreter','latex'); 




%--------------------- Plot C -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{1,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');

% Calibration
for i = 1:n_split
    % Calcula valor médio de SSE para primeira combinacao de dados
    Media = mean(EPR_out_CrossVal(i).Gen_Cal_r,2);
    % Obter intervalos de confiança (95%)
    [SSE.low, SSE.up] = confidence_intervals (EPR_out_CrossVal(i).Gen_Cal_r');
    
    j=i-0.15;
    plot([j-0.12, j+0.12],[SSE.up(1,DE.last),SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot([j-0.12, j+0.12],[SSE.low(1,DE.last),SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot([j,j],[Media(DE.last,1), SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot([j,j],[Media(DE.last,1), SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot(j,Media(DE.last,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_cal);hold on;
end

% Eval
for i = 1:n_split
    % Calcula valor médio de SSE para primeira combinacao de dados
    Media = mean(EPR_out_CrossVal(i).Gen_eval_r,2);
    % Obter intervalos de confiança (95%)
    [SSE.low, SSE.up] = confidence_intervals (EPR_out_CrossVal(i).Gen_eval_r');
    
    j=i+0.15;
    plot([j-0.12, j+0.12],[SSE.up(1,DE.last),SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot([j-0.12, j+0.12],[SSE.low(1,DE.last),SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot([j,j],[Media(DE.last,1), SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot([j,j],[Media(DE.last,1), SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot(j,Media(DE.last,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_eval);hold on;
end


set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
xlabel('Cross-validation batch','interpreter','latex','FontSize',fax);
ylabel('$r$','interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);


set(gca,'XLim',[1-0.5 i+0.5]);
xticks([1 2 3 4 5 6 7 8 9 10])
text(0.905,0.925,'(c)','units','normalized','FontSize',fax,'interpreter','latex'); 


%--------------------- Plot D -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{2,1},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');

% Calibration
for i = 1:n_split
    % Calcula valor médio de SSE para primeira combinacao de dados
    Media = mean(EPR_out_CrossVal(i).Gen_Cal_PBIAS,2);
    % Obter intervalos de confiança (95%)
    [SSE.low, SSE.up] = confidence_intervals (EPR_out_CrossVal(i).Gen_Cal_PBIAS');
    
    j=i-0.15;
    plot([j-0.12, j+0.12],[SSE.up(1,DE.last),SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot([j-0.12, j+0.12],[SSE.low(1,DE.last),SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot([j,j],[Media(DE.last,1), SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot([j,j],[Media(DE.last,1), SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_cal);hold on;
    plot(j,Media(DE.last,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_cal);hold on;
end

% Eval
for i = 1:n_split
    % Calcula valor médio de SSE para primeira combinacao de dados
    Media = mean(EPR_out_CrossVal(i).Gen_eval_PBIAS,2);
    % Obter intervalos de confiança (95%)
    [SSE.low, SSE.up] = confidence_intervals (EPR_out_CrossVal(i).Gen_eval_PBIAS');
    
    j=i+0.15;
    plot([j-0.12, j+0.12],[SSE.up(1,DE.last),SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot([j-0.12, j+0.12],[SSE.low(1,DE.last),SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot([j,j],[Media(DE.last,1), SSE.low(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot([j,j],[Media(DE.last,1), SSE.up(1,DE.last)],'k-','Linewidth',1.5,'Color',Color_eval);hold on;
    plot(j,Media(DE.last,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_eval);hold on;
end


set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1        );
xlabel('Cross-validation batch','interpreter','latex','FontSize',fax);
ylabel('PBIAS ($\%$)','interpreter','latex','FontSize',fax);
grid on
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',fax);


set(gca,'XLim',[1-0.5 i+0.5]);
xticks([1 2 3 4 5 6 7 8 9 10])
text(0.905,0.925,'(d)','units','normalized','FontSize',fax,'interpreter','latex'); 



%----------------------------- Sensitivity -------------------------------%
% compute independent variables
n_x = size(X_cal,2);

plotheight=12;%12
plotwidth=25;%8
subplotsx=n_x;
subplotsy=4;   
leftedge=1.5;
rightedge=0.5;   
topedge=0.5;
bottomedge=1.5;
spacex=1.8;
spacey=1.5;
fontsize=12;    
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

% setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'units','normalized','outerposition',[0.1 0 .7 .7])



Color_ini = [.7 .7 0];
Color_mid = [0 .7 .7];
Color_end = [.7 0 .7];

% Number of parameters
for i=1:n_x

%--------------------- Plot A -----------------------------%
% define position and initial settings for the axis
ax=axes('position',sub_pos{i,4},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
    
    
    % Number of MonteCarlo runs
    for j = 1: n_MonteCarlo 
        X_value(:,j) = X_sensitivity.X(:,i,i);
        Y_value(:,j) = EPR_out_CrossVal(1).Y_sens_MonteCarlo(:,i,j);
    end
    [Y_range.low, Y_range.up] = confidence_intervals (Y_value');
    Rang = Fill_Ranges(X_sensitivity.maxmin(:,i),Y_range.low,Y_range.up, [0.5 0.5 0.5]);hold on 
    P4 = plot(X_sensitivity.maxmin(:,i),mean(Y_value,2),'k--','LineStyle','-','Linewidth',2);  hold on
    
    alpha(Rang, 0.3);
    
    % Plot error bars
    X_min = X_sensitivity.maxmin(1,i);
    mid = round(size(X_cal,1)/2);
    X_mid = X_sensitivity.maxmin(mid,i);
    X_max = X_sensitivity.maxmin(end,i);
    Y_media = mean(Y_value,2);
    
    % Automize size of confidence inteerval
    if X_max-X_min <= 0.0001
        tam = 0.000001 ;
    elseif 0.0001< X_max-X_min && X_max-X_min < 1
        tam = 0.0075;
    else
        tam = 5;
    end
        
    plot([X_min-tam, X_min+tam],[Y_range.up(1,1),Y_range.up(1,1)],'k-','Linewidth',1.5,'Color',Color_ini);hold on;
    plot([X_min-tam, X_min+tam],[Y_range.low(1,1),Y_range.low(1,1)],'k-','Linewidth',1.5,'Color',Color_ini);hold on;
    plot([X_min,X_min],[Y_media(1,1), Y_range.low(1,1)],'k-','Linewidth',1.5,'Color',Color_ini);hold on;
    plot([X_min,X_min],[Y_media(1,1), Y_range.up(1,1)],'k-','Linewidth',1.5,'Color',Color_ini);hold on;
    plot(X_min,Y_media(1,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_ini);hold on;

    plot([X_mid-tam, X_mid+tam],[Y_range.up(1,mid),Y_range.up(1,mid)],'k-','Linewidth',1.5,'Color',Color_mid);hold on;
    plot([X_mid-tam, X_mid+tam],[Y_range.low(1,mid),Y_range.low(1,mid)],'k-','Linewidth',1.5,'Color',Color_mid);hold on;
    plot([X_mid,X_mid],[Y_media(mid,1), Y_range.low(1,mid)],'k-','Linewidth',1.5,'Color',Color_mid);hold on;
    plot([X_mid,X_mid],[Y_media(mid,1), Y_range.up(1,mid)],'k-','Linewidth',1.5,'Color',Color_mid);hold on;
    plot(X_mid,Y_media(mid,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_mid);hold on;
  
    plot([X_max-tam, X_max+tam],[Y_range.up(1,end),Y_range.up(1,end)],'k-','Linewidth',1.5,'Color',Color_end);hold on;
    plot([X_max-tam, X_max+tam],[Y_range.low(1,end),Y_range.low(1,end)],'k-','Linewidth',1.5,'Color',Color_end);hold on;
    plot([X_max,X_max],[Y_media(end,1), Y_range.low(1,end)],'k-','Linewidth',1.5,'Color',Color_end);hold on;
    plot([X_max,X_max],[Y_media(end,1), Y_range.up(1,end)],'k-','Linewidth',1.5,'Color',Color_end);hold on;
    plot(X_max,Y_media(end,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_end);hold on;

    

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

    if i==3
 
    plot([X_min-X_min*0.2, X_min+X_min*0.2],[Y_range.up(1,1),Y_range.up(1,1)],'k-','Linewidth',1.5,'Color',Color_ini);hold on;
    plot([X_min-X_min*0.2, X_min+X_min*0.2],[Y_range.low(1,1),Y_range.low(1,1)],'k-','Linewidth',1.5,'Color',Color_ini);hold on;
    plot([X_min,X_min],[Y_media(1,1), Y_range.low(1,1)],'k-','Linewidth',1.5,'Color',Color_ini);hold on;
    plot([X_min,X_min],[Y_media(1,1), Y_range.up(1,1)],'k-','Linewidth',1.5,'Color',Color_ini);hold on;
    plot(X_min,Y_media(1,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_ini);hold on;
   
        
    plot([X_mid-X_mid*0.2, X_mid+X_mid*0.2],[Y_range.up(1,mid),Y_range.up(1,mid)],'k-','Linewidth',1.5,'Color',Color_mid);hold on;
    plot([X_mid-X_mid*0.2, X_mid+X_mid*0.2],[Y_range.low(1,mid),Y_range.low(1,mid)],'k-','Linewidth',1.5,'Color',Color_mid);hold on;
    plot([X_mid,X_mid],[Y_media(mid,1), Y_range.low(1,mid)],'k-','Linewidth',1.5,'Color',Color_mid);hold on;
    plot([X_mid,X_mid],[Y_media(mid,1), Y_range.up(1,mid)],'k-','Linewidth',1.5,'Color',Color_mid);hold on;
    plot(X_mid,Y_media(mid,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_mid);hold on;
  
    plot([X_max-X_max*0.2, X_max+X_max*0.2],[Y_range.up(1,end),Y_range.up(1,end)],'k-','Linewidth',1.5,'Color',Color_end);hold on;
    plot([X_max-X_max*0.2, X_max+X_max*0.2],[Y_range.low(1,end),Y_range.low(1,end)],'k-','Linewidth',1.5,'Color',Color_end);hold on;
    plot([X_max,X_max],[Y_media(end,1), Y_range.low(1,end)],'k-','Linewidth',1.5,'Color',Color_end);hold on;
    plot([X_max,X_max],[Y_media(end,1), Y_range.up(1,end)],'k-','Linewidth',1.5,'Color',Color_end);hold on;
    plot(X_max,Y_media(end,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_end);hold on;
      
    set(gca,'XLim',[1e-10 1e-6]);
    set(gca,'XScale','log'); 
    
    end
    
    if i==1
        % Error bar for legend
            error1 = errorbar(-20,1,1,'o-','Color',Color_ini,...
                    'MarkerFaceColor','w','LineWidth',1.5,'LineStyle','none'); hold on;
            error2 = errorbar(-20,1,1,'o-','Color',Color_mid,...
                    'MarkerFaceColor','w','LineWidth',1.5,'LineStyle','none'); hold on;
            error3 = errorbar(-20,1,1,'o-','Color',Color_end,...
                    'MarkerFaceColor','w','LineWidth',1.5,'LineStyle','none'); hold on;
            leg1 = legend([Rang, P4, error1, error2, error3],['95$\%$ uncert.'], ...
                ['Mean'],['Min'],['Mid.'],['Max'],...
                    'location', 'southwest','Orientation','vertical');
            set(leg1,'Interpreter','latex');
            
            set(gca,'XLim',[0 0.3]);
    end
    
    
    
    
    %------------- Plot cross validation ----------------%
    %--------------------- Plot B -----------------------------%
    % define position and initial settings for the axis
    ax=axes('position',sub_pos{i,3},...
                'XGrid','off',...
                'FontSize',fontsize,...
                'Layer','top');  
            
    for c = 1:n_split
        for j = 1: n_MonteCarlo 
            X_value(:,j) = X_sensitivity.X(:,i,i);
            Y_value(:,j) = EPR_out_CrossVal(c).Y_sens_MonteCarlo(:,i,j);
        end
        [Y_range.low, Y_range.up] = confidence_intervals (Y_value');

        % organize dependente variable
       % X_min = X_sensitivity.maxmin(1,i);
       % mid = round(size(X_cal,1)/2);
       % X_mid = X_sensitivity.maxmin(mid,i);
       % X_max = X_sensitivity.maxmin(end,i);
        Y_media = mean(Y_value,2);

        % inserir outros X_med e Xmax
        plot([c-.2, c+.2],[Y_range.up(1,1),Y_range.up(1,1)],'k-','Linewidth',1.5,'Color',Color_ini);hold on;
        plot([c-.2, c+.2],[Y_range.low(1,1),Y_range.low(1,1)],'k-','Linewidth',1.5,'Color',Color_ini);hold on;
        plot([c,c],[Y_media(1,1), Y_range.low(1,1)],'k-','Linewidth',1.5,'Color',Color_ini);hold on;
        plot([c,c],[Y_media(1,1), Y_range.up(1,1)],'k-','Linewidth',1.5,'Color',Color_ini);hold on;
        plot(c,Y_media(1,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_ini);hold on;
    
    %    xlabel(cellstr(Names(i+1)),'interpreter','latex','FontSize',fax);
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
    
    set(gca,'XLim',[1-0.5 n_split+0.5]);
    xticks([1 2 3 4 5 6 7 8 9 10])

    end
    
    
    
    %--------------------- Plot C -----------------------------%
    % define position and initial settings for the axis
    ax=axes('position',sub_pos{i,2},...
                'XGrid','off',...
                'FontSize',fontsize,...
                'Layer','top');
    
    for c = 1:n_split
        for j = 1: n_MonteCarlo 
            X_value(:,j) = X_sensitivity.X(:,i,i);
            Y_value(:,j) = EPR_out_CrossVal(c).Y_sens_MonteCarlo(:,i,j);
        end
        [Y_range.low, Y_range.up] = confidence_intervals (Y_value');

        % organize dependente variable
       % X_min = X_sensitivity.maxmin(1,i);
        mid = round(size(X_cal,1)/2);
       % X_mid = X_sensitivity.maxmin(mid,i);
       % X_max = X_sensitivity.maxmin(end,i);
        Y_media = mean(Y_value,2);
        
        % inserir outros X_med e Xmax
        plot([c-.2, c+.2],[Y_range.up(1,mid),Y_range.up(1,mid)],'k-','Linewidth',1.5,'Color',Color_mid);hold on;
        plot([c-.2, c+.2],[Y_range.low(1,mid),Y_range.low(1,mid)],'k-','Linewidth',1.5,'Color',Color_mid);hold on;
        plot([c,c],[Y_media(mid,1), Y_range.low(1,mid)],'k-','Linewidth',1.5,'Color',Color_mid);hold on;
        plot([c,c],[Y_media(mid,1), Y_range.up(1,mid)],'k-','Linewidth',1.5,'Color',Color_mid);hold on;
        plot(c,Y_media(mid,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_mid);hold on;
    
    %    xlabel(cellstr(Names(i+1)),'interpreter','latex','FontSize',fax);
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
    
    set(gca,'XLim',[1-0.5 n_split+0.5]);
    xticks([1 2 3 4 5 6 7 8 9 10])

    end

    
    
     %--------------------- Plot D -----------------------------%
    % define position and initial settings for the axis
    ax=axes('position',sub_pos{i,1},...
                'XGrid','off',...
                'FontSize',fontsize,...
                'Layer','top');
    
    for c = 1:n_split

        for j = 1: n_MonteCarlo 
            X_value(:,j) = X_sensitivity.X(:,i,i);
            Y_value(:,j) = EPR_out_CrossVal(c).Y_sens_MonteCarlo(:,i,j);
        end
        [Y_range.low, Y_range.up] = confidence_intervals (Y_value');
        
        Y_media = mean(Y_value,2);
        
        % inserir outros X_med e Xmax
        plot([c-.2, c+.2],[Y_range.up(1,end),Y_range.up(1,end)],'k-','Linewidth',1.5,'Color',Color_end);hold on;
        plot([c-.2, c+.2],[Y_range.low(1,end),Y_range.low(1,end)],'k-','Linewidth',1.5,'Color',Color_end);hold on;
        plot([c,c],[Y_media(end,1), Y_range.low(1,end)],'k-','Linewidth',1.5,'Color',Color_end);hold on;
        plot([c,c],[Y_media(end,1), Y_range.up(1,end)],'k-','Linewidth',1.5,'Color',Color_end);hold on;
        plot(c,Y_media(end,1),'ko','MarkerFaceColor','w','MarkerSize',6,'Linewidth',1.5,'Color',Color_end);hold on;
    
    xlabel('Cross-validation batch','interpreter','latex','FontSize',fax);
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
    
    set(gca,'XLim',[1-0.5 n_split+0.5]);
    xticks([1 2 3 4 5 6 7 8 9 10])
    
    end

    
end
    











%------------------------ Sensitivity Adjustment ------------------------%
% compute independent variables


plotheight=12;%12
plotwidth=25;%8
subplotsx=n_x;
subplotsy=3;   
leftedge=1.5;
rightedge=0.5;   
topedge=0.5;
bottomedge=1.5;
spacex=1.8;
spacey=1.5;
fontsize=12;    
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

% setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'units','normalized','outerposition',[0.1 0 .7 .7])


% Number of parameters
for i=1:n_x

%--------------------- Plot A -----------------------------%
% define position and initial settings for the axis
%ax=axes('position',sub_pos{1,n_x+1-i},...
ax=axes('position',sub_pos{i,3},...
            'XGrid','off',...
            'FontSize',fontsize,...
            'Layer','top');
    
    
    % Number of MonteCarlo runs
    for j = 1: n_MonteCarlo 
        X_value(:,j) = X_sensitivity.X(:,i,i);
        Y_value(:,j) = EPR_out_CrossVal(1).Y_sens_MonteCarlo(:,i,j);
        P1 = plot(X_sensitivity.maxmin(:,i),EPR_out_CrossVal(1).Y_sens_MonteCarlo(:,i,j),'k--','Linewidth',0.5, 'Color', [.5 .5 .5]);  hold on
    end
    % Plot Selected model
    P2 = plot(X_sensitivity.maxmin(:,i),EPR_out_CrossVal(1).Y_sens_MonteCarlo(:,i,EPR_out_CrossVal(1).BEST),'k--','Linewidth',2, 'Color', [.0 .0 1]);  hold on
 
    % Plot Mean
    P3 = plot(X_sensitivity.maxmin(:,i),mean(Y_value,2),'k--','LineStyle','-','Linewidth',2);  hold on
    
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

    if i==2
        set(gca,'XScale','log'); 
    end
    
    if i==1
        leg1 = legend([P1, P2, P3],['MonteCarlo'],['Selection'],['Mean'],'location', 'southwest','Orientation','vertical');
             set(leg1,'Interpreter','latex');
    end
    
    %------------- Plot cross validation adjustment ----------------%
    %--------------------- Plot B -----------------------------%
    % define position and initial settings for the axis
    %ax=axes('position',sub_pos{2,n_x+1-i},...
    ax=axes('position',sub_pos{i,2},...
                'XGrid','off',...
                'FontSize',fontsize,...
                'Layer','top');  
            
    for c = 1:n_split
        
        
        % COMPUTE STATISTICS - DISTANCE FROM EACH MODEL AGAINST THE MEAN
        for j = 1: n_MonteCarlo 
            % Get MonteCarlo Values
            Y_value(:,j) = EPR_out_CrossVal(c).Y_sens_MonteCarlo(:,i,j);
        end
        
        % Get Mean
        Y_med = mean(Y_value,2);
        
        % Get statistics
        for j = 1: n_MonteCarlo
            [Out_stats] = statistic_model (Y_med, Y_value(:,j)) ;
            % Check negative values (only rsq, r and E_rel) and replace 
            if Out_stats.rsq < 0
                Out_stats.rsq = 0.01;
            end
            % Save statistic for all monte carlo models (j) and ind. variable
            model_stats(c).rsq(j,1) = Out_stats.rsq;
            model_stats(c).RMSE(j,1) = Out_stats.RMSE;
        end
        
    % Box plot    
    %h_sens = boxplot2([model_stats(c).rsq], c); hold on;   
    B1 = boxchart(ones(n_MonteCarlo,1)*c,model_stats(c).rsq, ...
                'MarkerSize',2,'BoxFaceColor','k','MarkerColor','k'); hold on;
    
                % Change settings
%     structfun(@(x) set(x, 'color', [0 0 0], ...
%         'markeredgecolor', [0 0 0]), h_sens);
%     set([h_sens.lwhis h_sens.uwhis], 'linestyle', '-','linewidth',1);
%     set([h_sens.box],'linewidth',1);
%     set([h_sens.med],'linewidth',1);
%     set(h_sens.out, 'marker', '.');
    
    Y_selected_rsq(1,c) = model_stats(c).rsq(EPR_out_CrossVal(c).BEST);
    Y_selected_RMSE(1,c) = model_stats(c).RMSE(EPR_out_CrossVal(c).BEST);   
    end
    
    % PLOTAR LINHA
    L1 = plot ([1:1:n_split],Y_selected_rsq, 'ko-.','color',[0 0 1], ...
        'MarkerFaceColor',[.8 .8 1],'Markersize',5,'LineWidth',1); hold on;
    
    ylabel('R$^2$','interpreter','latex','FontSize',fax);
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
    

    if i==1
        leg1 = legend([B1, L1],['MonteCarlo'],['Selection'],'location', 'southwest','Orientation','vertical');
             set(leg1,'Interpreter','latex');
    end
    

    set(gca,'XLim',[1-0.5 n_split+0.5]);
    xticks([1 2 3 4 5 6 7 8 9 10]);
    
    %--------------------- Plot C -----------------------------%
    % define position and initial settings for the axis
    %ax=axes('position',sub_pos{3,n_x+1-i},...
    ax=axes('position',sub_pos{i,1},...
                'XGrid','off',...
                'FontSize',fontsize,...
                'Layer','top');
    
    for c = 1:n_split
        
%         h_sens = boxplot2([model_stats(c).RMSE], c); hold on;    
% 
%             % Change settings
%         structfun(@(x) set(x, 'color', [0 0 0], ...
%             'markeredgecolor', [0 0 0]), h_sens);
%         set([h_sens.lwhis h_sens.uwhis], 'linestyle', '-','linewidth',1);
%         set([h_sens.box],'linewidth',1);
%         set([h_sens.med],'linewidth',1);
%         set(h_sens.out, 'marker', '.');


        boxchart(ones(n_MonteCarlo,1)*c,model_stats(c).RMSE, ...
            'MarkerSize',2,'BoxFaceColor','k','MarkerColor','k'); hold on;
    end
    
        % PLOTAR LINHA
    plot ([1:1:n_split],Y_selected_RMSE, 'ko-.','color',[0 0 1], ...
        'MarkerFaceColor',[.8 .8 1],'Markersize',5,'LineWidth',1); hold on;
    
    ylabel('RMSE','interpreter','latex','FontSize',fax);
    xlabel('Cross-validation batch','interpreter','latex','FontSize',fax);
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
    
    set(gca,'XLim',[1-0.5 n_split+0.5]);
    xticks([1 2 3 4 5 6 7 8 9 10]);

end    


































%----------------------- Training and testing data ----------------------%

% plot
f = figure;
f.Position = [100 100 900 300];
% Plot each cross-validation batch
for i = 1:n_split
    % Plot calibration
    h_cal = boxplot2([EPR_out_CrossVal(i).Y_cal], i-0.1); hold on;
    % Change settings
    structfun(@(x) set(x, 'color', Color_cal, ...
        'markeredgecolor', Color_cal), h_cal);
    set([h_cal.lwhis h_cal.uwhis], 'linestyle', '-','linewidth',1.5);
    set([h_cal.box],'linewidth',1.5);
    set([h_cal.med],'linewidth',4);
    set(h_cal.out, 'marker', '.');
    
    % Plot evaluation
    h_eval = boxplot2([EPR_out_CrossVal(i).Y_eval], i+0.1); hold on;
    % Change settings
    structfun(@(x) set(x, 'color', Color_eval, ...
        'markeredgecolor', Color_eval), h_eval);
    set([h_eval.lwhis h_eval.uwhis], 'linestyle', '-','linewidth',1.5);
    set([h_eval.box],'linewidth',1.5);
    set([h_eval.med],'linewidth',4);
    set(h_eval.out, 'marker', '.');
    
end

% Ranges
xlim([0.5 n_split+0.5])
ylim([0 42])
set(gca,'fontsize',fax)


%pl = line([0.5 5.5],[1 1],'Color','r','LineStyle','--');
%xticklabels({'0','30','120','330'})
box on
%xlim([Par_info.min(1) Par_info.max(1)])
xlabel('Cross-validation batch','interpreter','latex','FontSize',fax);
ylabel('OMC, ($\%$)','interpreter','latex','FontSize',fax);

% interpretacao em Latex para os eixos
set(gca,'TickLabelInterpreter','latex');

% legend
h_cal = rectangle('Position',[4,36.5,.4,2],...
                'Curvature',0,'EdgeColor',Color_cal,...
                'LineWidth',1.5,'LineStyle','-');hold on;
h_eval = rectangle('Position',[5.5, 36.5,.4,2],...
                'Curvature',0,'EdgeColor',Color_eval,...
                'LineWidth',1.5,'LineStyle','-');hold on;
h = rectangle('Position',[3.8,35,3,5],...
                'Curvature',0,'EdgeColor','k',...
                'LineWidth',.5,'LineStyle','-');hold on;
hline = line(NaN,NaN,'LineWidth',1,'LineStyle','-','Color',[1 0 0]);hold on;

text(4.5, 37,'Train','FontSize',fax,'interpreter','latex')
text(6, 37,'Test','FontSize',fax,'interpreter','latex')









%% Tables
delete SplitSampling.txt

fileID = fopen('SplitSampling.txt','a');
fprintf(fileID,'%12s\r\n','---------- Description of Cross-validation batches ----------');
fprintf(fileID,'%12s\r\n','EPR algorithm: CrosMonteCarlo MODEGA with optimum length');
  for i =1:n_split
    fprintf(fileID,'%12s\r\n','--------------------------------------------');
    fprintf(fileID,'%12s %6.4f\r\n','Cross-val Batch =',i);
    fprintf(fileID,'%12s\r\n','              Train                Test                Difference');
    fprintf(fileID,'%12s %6.4f %12s %6.4f %12s %6.4f\r\n','     Max =', max(EPR_out_CrossVal(i).Y_cal),'    ',  max(EPR_out_CrossVal(i).Y_eval),'    ',  ((max(EPR_out_CrossVal(i).Y_cal)-max(EPR_out_CrossVal(i).Y_eval))/max(EPR_out_CrossVal(i).Y_cal))*100);
    fprintf(fileID,'%12s %6.4f %12s %6.4f %12s %6.4f\r\n','     Min =', min(EPR_out_CrossVal(i).Y_cal),'       ',  min(EPR_out_CrossVal(i).Y_eval),'    ',  ((min(EPR_out_CrossVal(i).Y_cal)-min(EPR_out_CrossVal(i).Y_eval))/min(EPR_out_CrossVal(i).Y_cal))*100);
    fprintf(fileID,'%12s %6.4f %12s %6.4f %12s %6.4f\r\n','    Mean =', mean(EPR_out_CrossVal(i).Y_cal),'    ',  mean(EPR_out_CrossVal(i).Y_eval),'    ',  ((mean(EPR_out_CrossVal(i).Y_cal)-mean(EPR_out_CrossVal(i).Y_eval))/mean(EPR_out_CrossVal(i).Y_cal))*100);
    fprintf(fileID,'%12s %6.4f %12s %6.4f %12s %6.4f\r\n','Std. dev.=', std(EPR_out_CrossVal(i).Y_cal),'       ',  std(EPR_out_CrossVal(i).Y_eval),'    ',  ((std(EPR_out_CrossVal(i).Y_cal)-std(EPR_out_CrossVal(i).Y_eval))/std(EPR_out_CrossVal(i).Y_cal))*100);
  end
fclose(fileID);


delete MathematicalStructures.txt

fileID = fopen('MathematicalStructures.txt','a');
fprintf(fileID,'%12s\r\n','---------- Display the selected equations of the Cross-validation batches ----------');
fprintf(fileID,'%12s\r\n','EPR algorithm: CrosMonteCarlo MODEGA with optimum length');
fprintf(fileID,'%12s\r\n',' Pop = Exponent values');
fprintf(fileID,'%12s\r\n',' Par = bias + terms values');
fprintf(fileID,'%12s %6.0f\r\n','Number of batches =         ',n_split);
fprintf(fileID,'%12s %6.0f\r\n','Number of terms (a0 + aj) = ',ndega+1);
fprintf(fileID,'%12s %6.0f\r\n','Number variables (X) =      ',size(X,2));
  for i =1:n_split  % number of split samplings
    fprintf(fileID,'%12s\r\n','--------------------------------------------');
    fprintf(fileID,'%12s %6.0f\r\n','Cross-val Batch =',i);
    fprintf(fileID,'%12s','Pop =');
    for j = 1:ndega % number of terms
        fprintf(fileID,'%6.2f', EPR_out_CrossVal(i).bestpop(j,:));
    end
    fprintf(fileID,'%12s\r\n','              ');
    fprintf(fileID,'%12s %6.2f','Par = ', EPR_out_CrossVal(i).best_par(1,1) );
    for j = 2:ndega+1 % number of terms
        fprintf(fileID,'%6.2f', EPR_out_CrossVal(i).best_par(1,j));
    end    
    fprintf(fileID,'%12s\r\n','              ');    
    %fprintf(fileID,'%12s\r\n','              Train                Test                Difference');
    %fprintf(fileID,'%12s %6.4f %12s %6.4f %12s %6.4f\r\n','     Max =', max(EPR_out_CrossVal(i).Y_cal),'    ',  max(EPR_out_CrossVal(i).Y_eval),'    ',  ((max(EPR_out_CrossVal(i).Y_cal)-max(EPR_out_CrossVal(i).Y_eval))/max(EPR_out_CrossVal(i).Y_cal))*100);
    %fprintf(fileID,'%12s %6.4f %12s %6.4f %12s %6.4f\r\n','     Min =', min(EPR_out_CrossVal(i).Y_cal),'       ',  min(EPR_out_CrossVal(i).Y_eval),'    ',  ((min(EPR_out_CrossVal(i).Y_cal)-min(EPR_out_CrossVal(i).Y_eval))/min(EPR_out_CrossVal(i).Y_cal))*100);
    %fprintf(fileID,'%12s %6.4f %12s %6.4f %12s %6.4f\r\n','    Mean =', mean(EPR_out_CrossVal(i).Y_cal),'    ',  mean(EPR_out_CrossVal(i).Y_eval),'    ',  ((mean(EPR_out_CrossVal(i).Y_cal)-mean(EPR_out_CrossVal(i).Y_eval))/mean(EPR_out_CrossVal(i).Y_cal))*100);
    %fprintf(fileID,'%12s %6.4f %12s %6.4f %12s %6.4f\r\n','Std. dev.=', std(EPR_out_CrossVal(i).Y_cal),'       ',  std(EPR_out_CrossVal(i).Y_eval),'    ',  ((std(EPR_out_CrossVal(i).Y_cal)-std(EPR_out_CrossVal(i).Y_eval))/std(EPR_out_CrossVal(i).Y_cal))*100);
  end
fclose(fileID);

                
    %---------------------  Plot uncertainty and average  --------------------%
    %[p]= plot_MonteCarlo_uncertainty(ndega, DE.last, n_MonteCarlo, EPR_out_MODEGA,EPR_out_MonteCarlo_MODEGA,EPR_eval_MonteCarlo_MODEGA);

    %----------  Plot sensitivity uncertainty, mean and best model -----------%

    % Plot sensitivity and output mean values of sensitivity
    %[Y_mean_sensitivity] = plot_sensitivity_ranges(n_MonteCarlo,X_sensitivity, EPR_out_MonteCarlo_MODEGA, X_cal, Y_cal, Names);

    % Function to select best model and plot 
    %[BEST]=model_selection(n_MonteCarlo,X_sensitivity,EPR_out_MonteCarlo_MODEGA, X_cal, Y_cal, Names);
