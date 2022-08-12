%--------------------------------------------------------------------------
% Learnin Koopman Operator
% Ref: Representer Theorem for Learning Koopman Operators 
% Link: https://arxiv.org/abs/2208.01681
%
% This code illustrates the results of evaluating the impact of different 
% regularization terms and constraint in learning problem 
%       min_K   E(K) + lambda * R(K),
%       s.t.    K in C,
% where E(.) is the sum squared error loss function, C is a given constraint, 
% R(.) is the regularization term. For C and R, we have the following cases:
%   1. R: square of operator norm of K, i.e.,  R(K) = ||K||^2 
%   2. R: square of Frobenius norm of K, i.e.,  R(K) = ||K||_F^2 
%   3. R: nuclear norm of K, i.e.,  R(K) = ||K||_* 
%   4. C: rank of K, i.e.,  C = {K | rank(K) <= r} 
%   5. R&C: if R=0 and C = L(K) (the whole space of bounded operators), the
%      program is equivalent to the EDMD approach
% For generating the data, two trajectories of the Van der Pol system is 
% considered. More details are provided in Example 2 in the reference above.  
%
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
close all;  clear all;  clc;
%--------------------------------------------------------------------------
disp('-------------------------------------------------------------------')
disp('Learnin Koopman Operator:                                          ')
disp('Ref: Representer Theorem for Learning Koopman Operators            ') 
disp('Link: https://arxiv.org/abs/2208.01681                             ')
disp('                                                                   ')
disp('This code illustrates the comparison results of employing different') 
disp('regularization terms and constraint in the learning problem        ') 
disp('       min_K   E(K) + lambda * R(K),                               ')
disp('       s.t.    K in C,                                             ')
disp('where E(.) is the sum squared error loss function, C corresponds to')
disp('given constraint, R(.) is the regularization term. For C and R, we ')
disp('consider the following cases:                                      ')
disp('   1. R: square of operator norm of K, i.e., R(K) = ||K||^2        ') 
disp('   2. R: square of Frobenius norm of K, i.e., R(K) = ||K||_F^2     ') 
disp('   3. R: nuclear norm of K, i.e.,  R(K) = ||K||_*                  ') 
disp('   4. C: rank of K, i.e.,  C = {K | rank(K) <= r}                  ') 
disp('   5. R&C: if R = 0 and C = L(K) (the space of bounded operators), ') 
disp('      the program is equivalent to the EDMD approach.              ')
disp('For generating the data, two trajectories of the Van der Pol system')
disp('is considered. More details and discussion are provided in Example ')
disp('2 in the reference above.                                          ')
disp('                                                                   ')
disp('Mohammad Khosravi                                                  ')
disp('Email: mohammad.khosravi@tudelft.nl                                ')
disp('Delft Center for Systems and Control (DCSC)                        ')
disp('Delft University of Technology (TU Delft)                          ') 
disp('August 2022                                                        ')
disp('-------------------------------------------------------------------')
disp(' ')

LW = 'linewidth';       FS = 'FontSize';        MS = 'MarkerSize';
LOC = 'Location';       INT = 'Interpreter';    LX = 'Latex';   
DN = 'DisplayName';     
OLS = 'OutlierSize';    FG = 'FactorGap';

CL = 'Color';
blue_CL     = [0 0.45 0.74];    sky_CL      = [0 0.75 0.75];
olive_CL    = [0.75 0.75 0];    orange_CL   = [0.85 0.33 0.1];
purple_CL   = [0.75 0 0.75];

%--------------------------------------------------------------------------
% loading data
N = 120;
err_MSE_edmd_10 = zeros(N,1);
err_MSE_edmd_20 = zeros(N,1);
err_MSE_edmd_30 = zeros(N,1);
err_MSE_fro_10  = zeros(N,1);
err_MSE_fro_20  = zeros(N,1);
err_MSE_fro_30  = zeros(N,1);
err_MSE_opr_10  = zeros(N,1);
err_MSE_opr_20  = zeros(N,1);
err_MSE_opr_30  = zeros(N,1);
err_MSE_nuc_10  = zeros(N,1);
err_MSE_nuc_20  = zeros(N,1);
err_MSE_nuc_30  = zeros(N,1);
err_MSE_rnk_10  = zeros(N,1);
err_MSE_rnk_20  = zeros(N,1);
err_MSE_rnk_30  = zeros(N,1);

dBvalx = 10;
for i = 1:N
    MC_num = i-1;
    MC_num_str = num2str(MC_num + 1000); MC_num_str = MC_num_str(2:end);
    fn = ['data_',MC_num_str,'_',num2str(dBvalx),'.mat'];
    load(fn)
    
    err_MSE_edmd_10(i) = err_MSE_edmd;
    err_MSE_fro_10(i)  = err_MSE_fro;
    err_MSE_opr_10(i)  = err_MSE_opr;
    err_MSE_nuc_10(i)  = err_MSE_nuc;
    err_MSE_rnk_10(i)  = err_MSE_rnk;
end

dBvalx = 20;
for i = 1:N
    MC_num = i-1;
    MC_num_str = num2str(MC_num + 1000); MC_num_str = MC_num_str(2:end);
    fn = ['data_',MC_num_str,'_',num2str(dBvalx),'.mat'];
    load(fn)
    
    err_MSE_edmd_20(i) = err_MSE_edmd;
    err_MSE_fro_20(i)  = err_MSE_fro;
    err_MSE_opr_20(i)  = err_MSE_opr;
    err_MSE_nuc_20(i)  = err_MSE_nuc;
    err_MSE_rnk_20(i)  = err_MSE_rnk;
end

dBvalx = 30;
for i = 1:N
    MC_num = i-1;
    MC_num_str = num2str(MC_num + 1000); MC_num_str = MC_num_str(2:end);
    fn = ['data_',MC_num_str,'_',num2str(dBvalx),'.mat'];
    load(fn)
    
    err_MSE_edmd_30(i) = err_MSE_edmd;
    err_MSE_fro_30(i)  = err_MSE_fro;
    err_MSE_opr_30(i)  = err_MSE_opr;
    err_MSE_nuc_30(i)  = err_MSE_nuc;
    err_MSE_rnk_30(i)  = err_MSE_rnk;
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% plotting the boxplots of the results

% Font sizes:
tFS   = 11;         % title font sizes
ylFS  = 11;         % ylable font sizes
xlFS  = 11;         % xlable font sizes
lgFS  = 12;         % legend font sizes
xtlFS = 10;         % x-tick font sizes
ytlFS = 10;         % y-tick font sizes
ols_value = 3;      % boxplot: linewidth for outliers +
Width_value = 0.5;  % boxplot: width

f = figure();
f.Position = [450 450 1200 300]; 

    %----------------------------------------------------------------------
    
    subplot(1,5,2)
    boxplot([err_MSE_opr_10,err_MSE_opr_20,err_MSE_opr_30],...
        OLS,ols_value,'Widths',Width_value)
    
    grid on
    ax = gca;
    ax.YScale ='log';
    ax.YLim = [1,2e6];
    ax.TickLabelInterpreter = LX;
    ax.XTickLabel(1) = {'10dB'};
    ax.XTickLabel(2) = {'20dB'};
    ax.XTickLabel(3) = {'30dB'};
    ax.YTick = [1 1e1 1e2 1e3 1e4 1e5 1e6];
    ax.YTickLabel(2) = {''};
    ax.YTickLabel(4) = {''};
    ax.YTickLabel(6) = {''};
    xtl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',xtl,'fontsize',xtlFS,'FontName','Times New Roman')
    t = title('$\|\cdot\|$ reg.', INT,LX);
    t.FontSize = tFS;
    hsp2 = get(gca, 'Position');
    
    %----------------------------------------------------------------------
    
    subplot(1,5,1)
    boxplot([err_MSE_edmd_10,err_MSE_edmd_20,err_MSE_edmd_30],...
        OLS,ols_value,'Widths',Width_value)
    
    grid on
    ax = gca;
    ax.YScale ='log';
    ax.YLim = [1,2e6];
    ax.TickLabelInterpreter = LX;
    ax.XTickLabel(1) = {'10dB'};
    ax.XTickLabel(2) = {'20dB'};
    ax.XTickLabel(3) = {'30dB'};
    ax.YTick = [1 1e1 1e2 1e3 1e4 1e5 1e6];
    ax.YTickLabel(2) = {''};
    ax.YTickLabel(4) = {''};
    ax.YTickLabel(6) = {''};
    xtl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',xtl,'fontsize',xtlFS,'FontName','Times New Roman')
    t = title('EDMD', INT,LX);
    t.FontSize = tFS;   
    yl = ylabel('mean squared error', INT,LX);
    yl.FontSize = ylFS;
    hsp1 = get(gca, 'Position') - [0.01 0  0 0];                
    set(gca, 'Position', [hsp1(1:2)  hsp2(3:4)]);
    
    %----------------------------------------------------------------------
    
    subplot(1,5,3)
    boxplot([err_MSE_fro_10,err_MSE_fro_20,err_MSE_fro_30],...
        OLS,ols_value,'Widths',Width_value)
    
    grid on
    ax = gca;
    ax.YScale ='log';
    ax.YLim = [1,2e6];
    ax.TickLabelInterpreter = LX;
    ax.XTickLabel(1) = {'10dB'};
    ax.XTickLabel(2) = {'20dB'};
    ax.XTickLabel(3) = {'30dB'};
    ax.YTick = [1 1e1 1e2 1e3 1e4 1e5 1e6];
    ax.YTickLabel(2) = {''};
    ax.YTickLabel(4) = {''};
    ax.YTickLabel(6) = {''};
    xtl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',xtl,'fontsize',xtlFS,'FontName','Times New Roman')
    t = title('$\|\cdot\|_{\mathrm{F}}$ reg.', INT,LX);
    t.FontSize = tFS;
    hsp3 = get(gca, 'Position');                
    set(gca, 'Position', [hsp3(1:2)  hsp2(3:4)]);
    
    %----------------------------------------------------------------------
    
    subplot(1,5,4)
    grid on
    boxplot([err_MSE_nuc_10,err_MSE_nuc_20,err_MSE_nuc_30],...
        OLS,ols_value,'Widths',Width_value)
    
    grid on
    ax = gca;
    ax.YScale ='log';
    ax.YLim = [1,2e6];
    ax.TickLabelInterpreter = LX;
    ax.XTickLabel(1) = {'10dB'};
    ax.XTickLabel(2) = {'20dB'};
    ax.XTickLabel(3) = {'30dB'};
    ax.YTick = [1 1e1 1e2 1e3 1e4 1e5 1e6];
    ax.YTickLabel(2) = {''};
    ax.YTickLabel(4) = {''};
    ax.YTickLabel(6) = {''};
    xtl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',xtl,'fontsize',xtlFS,'FontName','Times New Roman')
    t = title('$\|\cdot\|_{*}$ reg.', INT,LX);
    t.FontSize = tFS;
    hsp4 = get(gca, 'Position');                
    set(gca, 'Position', [hsp4(1:2)  hsp2(3:4)]);
    
    %----------------------------------------------------------------------
    
    subplot(1,5,5)
    boxplot([err_MSE_rnk_10,err_MSE_rnk_20,err_MSE_rnk_30],...
        OLS,ols_value,'Widths',Width_value)
    
    grid on
    ax = gca;
    ax.YScale ='log';
    ax.YLim = [1,2e6];
    ax.TickLabelInterpreter = LX;
    ax.XTickLabel(1) = {'10dB'};
    ax.XTickLabel(2) = {'20dB'};
    ax.XTickLabel(3) = {'30dB'};
    ax.YTick = [1 1e1 1e2 1e3 1e4 1e5 1e6];
    ax.YTickLabel(2) = {''};
    ax.YTickLabel(4) = {''};
    ax.YTickLabel(6) = {''};
    xtl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',xtl,'fontsize',xtlFS,'FontName','Times New Roman')
    t = title('rank const.', INT,LX);
    t.FontSize = tFS;
    hsp5 = get(gca, 'Position');                
    set(gca, 'Position', [hsp5(1:2)  hsp2(3:4)]);
