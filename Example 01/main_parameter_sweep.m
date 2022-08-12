%--------------------------------------------------------------------------
% Learnin Koopman Operator
% Ref: Representer Theorem for Learning Koopman Operators 
% Link: https://arxiv.org/abs/2208.01681
%
% Example: regularization weight sweeping (Example 1 in ref)
% In this example, we consider the following learning problem 
%   min_K E(K) + lambda * R(K),
% where E(.) is the sum squared error loss function and R(.) is one of the 
% following regularization terms:
%   1. square of operator norm of K, i.e.,  R(K) = ||K||^2 
%   2. square of Frobenius norm of K, i.e.,  R(K) = ||K||_F^2 
%   3. nuclear norm of K, i.e.,  R(K) = ||K||_* 
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
disp('This code concerns sweeping the regularization weight lambda in the') 
disp('learning Koopman operator problem introduced in the reference above')  
disp('(Example 1). We consider the following learning problem            ')
disp('      min_K E(K) + lambda * R(K),                                  ')
disp('where E(.) is the sum squared error loss function and, R(.) is one ')
disp('of the following regularization terms:                             ')
disp('   1. square of operator norm of K, i.e., R(K) = ||K||^2,          ')
disp('   2. square of Frobenius norm of K, i.e., R(K) = ||K||_F^2,       ')
disp('   3. nuclear norm of K, i.e., R(K) = ||K||_*.                     ')
disp('We sweep lambda from 0 to infinite, and solve the learning problem.')
disp('                                                                   ')
disp('Mohammad Khosravi                                                  ')
disp('Email: mohammad.khosravi@tudelft.nl                                ')
disp('Delft Center for Systems and Control (DCSC)                        ')
disp('Delft University of Technology (TU Delft)                          ')
disp('August 2022                                                        ')
disp('-------------------------------------------------------------------')
disp('The program is started!')
format shortg
c = clock;
disp(['Start time: ',num2str(c(1)),'-',num2str(c(2))','-',num2str(c(3)),' ',...
    num2str(c(4)),':',num2str(c(5))])
disp(' ')
%--------------------------------------------------------------------------

LW = 'linewidth';       FS = 'FontSize';        MS = 'MarkerSize';
LOC = 'Location';       INT = 'Interpreter';    LX = 'Latex';   
DN = 'DisplayName';     
OLS = 'OutlierSize';    FG = 'FactorGap';

CL = 'Color';
blue_CL     = [0 0.45 0.74];    sky_CL      = [0 0.75 0.75];
olive_CL    = [0.75 0.75 0];    orange_CL   = [0.85 0.33 0.1];
purple_CL   = [0.75 0 0.75];

%--------------------------------------------------------------------------
% adding the path of support files
currentFolder = pwd;
addpath([currentFolder(1:end-10),'support_codes'])

%--------------------------------------------------------------------------
% Generating tajectory of the system
% range:
x1_min = -1;     x1_max = 1;
x2_min = -0.5;   x2_max = 1.5;

% initial point:
Xi0 = [1, 0]; 

% number of trajectories
nt = size(Xi0,1);

% time range
T  = 0:1:60;

% trajectory data
X  = [];   % X:   [x(0) x(1) x(2) ... x(nS-1)      ] (for single trajectory)
Xp = [];   % Xp:  [     x(1) x(2) ... x(nS-1) x(nS)] (for single trajectory)
XXp = [];  % XXp: [x(0) x(1) x(2) ... x(nS-1) x(nS)] (for single trajectory)

% some index set for the trjaectories
idx_X = [];     idx_Xp = [];    idx_X_01 = [];  idx_Xp_01 = []; idx_XXp = [];


for i = 1:nt
    % initiating dynamical system
    Xi = Xi0(i,:);   

    idx_X_j = [];
    idx_Xp_j = [];

    % generating the system i-th trajectory
        for j=1:length(T)-1
            % next sample
            Xi(j+1,:) = (DiffE_2D(Xi(j,:)))';
            
            % In case we want to avoid similar trajectory point, we can set         
            % eps_X to small value and use following lines:
            % eps_X = norm(Xi(end-1,:)-Xi(end,:));
            % if eps_X < eps_X_lim
            %     break;
            % end
            
            % collecting the data point
            idx_X_j  = [idx_X_j;  i j-1];
            idx_Xp_j = [idx_Xp_j; i j];
        end
    
    % concatenate data points of trajectories    
    X   = [X;  Xi(1:end-1,:)];
    Xp  = [Xp; Xi(2:end,:)];
    XXp = [XXp; Xi(1:end,:)];
    
    % updating index sets
    idx_X   = [idx_X;    idx_X_j];
    idx_Xp  = [idx_Xp;   idx_Xp_j];
    idx_XXp = [idx_XXp;  idx_X_j; i j];    
    idx_X_01   = [idx_X_01;   ones(j,1); 0];
    idx_Xp_01  = [idx_Xp_01;  0;  ones(j,1)];
end

% We have:
% XXp(idx_X_loc,:)  = X  or X  = XXp(idx_X_loc,:) 
% XXp(idx_Xp_loc,:) = Xp or Xp = XXp(idx_Xp_loc,:) 
idx_X_loc  = find(idx_X_01  == 1);
idx_Xp_loc = find(idx_Xp_01 == 1);

% to save data points:
% save('data_DiffE_parameter_sweep.mat',...
%     'Xi0','nt', ...
%     'X','Xp','XXp',...
%     'idx_X','idx_Xp',...
%     'idx_X_01','idx_Xp_01','idx_XXp',...
%     'idx_X_loc','idx_Xp_loc')

%--------------------------------------------------------------------------
% Choosing kernel type and parameters
kernel_type = 'SE';
theta_k = [1 1];

%--------------------------------------------------------------------------
% The observales are defined as 
%       g_l(.) = k(p_l, .),     l = 1,..., nG,
% where k is the kernel function. The matrix P collects p_1,...,p_nG, i.e.,
% P = [p_1',p_2',...,p_nG']'
P = [    1     0;
         1     1;   
         0     1;];

% KG is the Gram matrix of p_1,...,p_nG, i.e., KG_(i,j) = k(p_i,p_j)    
KG  = kernel_fun(  P, P, theta_k, kernel_type);

% making sure that G is positive definite and generating its square root Gr.
eps_KG = 1e-4;
[KG,KGr] = mxPD(KG,eps_KG);
nG = size(KG,1);

% % trajectory data
% X  = XXp(idx_X_loc,:);
% Xp = XXp(idx_Xp_loc,:);

% Calculating matrix KXXpP with entries KXXpP_ij = k(x(i),p_j), forci = 0, 
% ... , nS and  j = 1, ... , nG
KXXpP = kernel_fun(XXp, P, theta_k, kernel_type);

% According to definition of Y, we can take it as a sub-matrix of KXXpP
Y = KXXpP(idx_Xp_loc,:);

% Calculating matrix KPG with entries KPG_ij = k(x(i),p_j), for i = 0, ... 
% , nS and  j = 1, ... , nG. Note that KPG is a sub-matrix of KXXpP
KPG = kernel_fun(X,P,theta_k,kernel_type);

% Here, we calculate matrices PGinv_GV = inv(G)*KPG and GV = KPG*inv(G)*KPG. 
% For more details, please see Remark 4 in the reference.
PGinv_GV = (eye(size(KG))/KG) * (KPG');
GV = KPG * PGinv_GV;

% To unify the notations, we use Z.
Z = GV; 

% making sure that Z is positive definite and generating its square root Zr.
eps_Z = 1e-2;
[Z,Zr] = mxPD(Z,eps_Z);
nZ = size(Z,1);

%--------------------------------------------------------------------------
% parameter sweep loop ----------------------------------------------------

% If you have installed cvx and Mosek, change the following variables to 1,
% this will provide more exact solutions.
cvx_solve_opr = 0;
cvx_solve_nuc = 0;

% initializing the loop
n_lambda = 250;
eig_fro = zeros(nG,n_lambda);
eig_opr = zeros(nG,n_lambda);
eig_nuc = zeros(nG,n_lambda);
eig_rnk = zeros(nG,n_lambda);

lambda_root = 1.2;
for i = 1:n_lambda
lambda = lambda_root^(-0.5 * (n_lambda-1)+i);
%--------------------------------------------------------------------------
% Operator norm regularization:
% We use change-of-variable B = Zr * A * Gr, solve the following CONVEX
% program
%       min_K ||Zr * B * Gr - Y||_F^2 + lambda * ||B||^2,
% and set A = inv(Zr) * B * inv(Gr). Note that without cvx + Mosek, we have
% some partial numerical inexactness in the solutions.

if cvx_solve_opr
    cvx_solver mosek
    cvx_begin quiet
        variable B(nZ,nG)
        variable opr_bnd(1)
        minimize(sum_square(vec(Zr * B * KGr - Y)) +  lambda * opr_bnd^2)
        subject to
        norm(B) <= opr_bnd
    cvx_end
    A_star_opr = (eye(nZ)/Zr) * B * (eye(nG)/KGr);
    
else
    [AB_star,~,~,~,~,~,~] = opt_sol_opr(Zr,KGr,Y,lambda);
    A_star_opr = AB_star;
end

C_star_opr = PGinv_GV * A_star_opr;
eig_opr(:,i) = eig(C_star_opr*KG);
%--------------------------------------------------------------------------
% Frobenius norm regularization:
% We would like to solve the following CONVEX programm
%       min_K ||Z * B * G - Y||_F^2 + lambda * ||Zr * B * Gr||_F^2,
% which has a closed form solution.

A_star_fro = opt_sol_fro(Z,KG,Y,lambda);
C_star_fro = PGinv_GV * A_star_fro;
eig_fro(:,i) = eig(C_star_fro*KG);

%--------------------------------------------------------------------------
% Nuclear norm regularization:
% We use change-of-variable B = Zr * A * Gr, solve the following CONVEX
% program
%       min_K ||Zr * B * Gr||_F^2 + lambda * ||B||_*,
% and set A = inv(Zr) * B * inv(Gr). Note that without cvx + Mosek, we have
% some partial numerical inexactness in the solutions.

if cvx_solve_nuc
    cvx_solver mosek
    cvx_begin quiet
        variable B(nZ,nG)
        minimize(sum_square(vec(Zr * B * KGr - Y)) + lambda * norm_nuc(B))
    cvx_end
    B_cvx = B;
    A_star_nuc = (eye(nZ)/Zr) * B * (eye(nG)/KGr);
else
    [AB_star,~,~,~,~,~,~] = opt_sol_nuc(Zr,KGr,Y,lambda);
    A_star_nuc = AB_star;
end
C_star_nuc = PGinv_GV * A_star_nuc;
eig_nuc(:,i) = eig(C_star_nuc*KG);

%--------------------------------------------------------------------------
perc = num2str(2000-floor(i/n_lambda*1000));
perc =[perc(2:3),'.',perc(4)];
disp([perc,'% of the parameter sweep is remained!'])

end
%--------------------------------------------------------------------------
% Extended dynamic mode decomposition (EDMD):
% We would like to solve the following CONVEX programm
%       min_K ||Z * B * G - Y||_F^2,
% which is the case of Frobenius norm with lambda = 0.

A_star_fro = opt_sol_fro(Z,KG,Y,0);
C_star_fro = PGinv_GV * A_star_fro;
eig_EDMD(:,i) = eig(C_star_fro*KG);

abs_eig_fro = sort(abs(eig_fro))';
abs_eig_opr = sort(abs(eig_opr))';
abs_eig_nuc = sort(abs(eig_nuc))';
abs_eig_rnk = sort(abs(eig_rnk))';
%--------------------------------------------------------------------------
% plotting the results

% Font sizes:
tFS   = 13;    % title font sizes
ylFS  = 11;    % ylable font sizes
xlFS  = 13;    % xlable font sizes
lgFS  = 12;    % legend font sizes
xtlFS = 10;    % x-tick font sizes

Lambda = lambda_root.^(-0.5 * (n_lambda-1)+(1:n_lambda));

f = figure();
f.Position = [100 100  700 570]; 

    subplot(3,1,1)
    hold on; box on; axis tight; grid on; 
    plot(Lambda(1)*ones(3,1),abs(eig_EDMD),'k o',LW,2.5,MS,3.5)
    plot(Lambda,abs_eig_opr(:,1),':',LW,1.5,CL, blue_CL)
    plot(Lambda,abs_eig_opr(:,2),':',LW,1.5,CL, orange_CL)
    plot(Lambda,abs_eig_opr(:,3),':',LW,1.5,CL, sky_CL)

    ax = gca;
    ax.YLim = [0,1.05];
    ax.TickLabelInterpreter = LX;
    ax.XScale ='log';     
    xtl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',xtl,'fontsize',xtlFS,'FontName','Times New Roman')
    ax.YTick = [0 0.25 0.5 0.75 1];
    ax.YTickLabel(2) = {''};
    ax.YTickLabel(4) = {''};
    t = title('Regularization: $\|\cdot\|^2$', INT,LX);
    t.FontSize = tFS;
    lg = legend('EDMD','$\gamma_1$','$\gamma_2$','$\gamma_3$',INT,LX);
    lg.FontSize = lgFS;
    hsp1 = get(gca, 'Position');
    %----------------------------------------------------------------------
    subplot(3,1,2)
    hold on; box on; axis tight; grid on;
    plot(Lambda(1)*ones(3,1),abs(eig_EDMD),'k o',LW,2.5,MS,3.5)
    plot(Lambda,abs_eig_fro(:,1),':',LW,1.5,CL, blue_CL)
    plot(Lambda,abs_eig_fro(:,2),':',LW,1.5,CL, orange_CL)
    plot(Lambda,abs_eig_fro(:,3),':',LW,1.5,CL, sky_CL)

    ax = gca;
    ax.YLim = [0,1.05];
    ax.TickLabelInterpreter = LX;
    ax.XScale ='log';
    xtl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',xtl,'fontsize',xtlFS,'FontName','Times New Roman')
    ax.YTick = [0 0.25 0.5 0.75 1];
    ax.YTickLabel(2) = {''};
    ax.YTickLabel(4) = {''};
    t = title('Regularization: $\|\cdot\|_{\mathrm{F}}^2$', INT,LX);
    t.FontSize = tFS;
    lg = legend('EDMD','$\gamma_1$','$\gamma_2$','$\gamma_3$',INT,LX);
    lg.FontSize = lgFS;
    hsp2 = get(gca, 'Position') + [0 0.01  0 0];                
    set(gca, 'Position', [hsp2(1:2)  hsp1(3:4)]);   
    %----------------------------------------------------------------------
    subplot(3,1,3)
    hold on; box on; axis tight; grid on;
    plot(Lambda(1)*ones(3,1),abs(eig_EDMD),'k o',LW,2.5,MS,3.5)
    plot(Lambda,abs_eig_nuc(:,1),':',LW,1.5,CL, blue_CL)
    plot(Lambda,abs_eig_nuc(:,2),':',LW,1.5,CL, orange_CL)
    plot(Lambda,abs_eig_nuc(:,3),':',LW,1.5,CL, sky_CL)
    
    ax = gca;
    ax.YLim = [0,1.05];
    ax.TickLabelInterpreter = LX;
    ax.XScale ='log'; 
    xtl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',xtl,'fontsize',xtlFS,'FontName','Times New Roman')
    ax.YTick = [0 0.25 0.5 0.75 1];
    ax.YTickLabel(2) = {''};
    ax.YTickLabel(4) = {''};
    t = title('Regularization: $\|\cdot\|_{*}$', INT,LX);
    t.FontSize = tFS;
    lg = legend('EDMD','$\gamma_1$','$\gamma_2$','$\gamma_3$',INT,LX);
    lg.FontSize = lgFS;
    hsp3 = get(gca, 'Position') + [0 0.02  0 0];                
    set(gca, 'Position', [hsp3(1:2)  hsp1(3:4)]);
    xl = xlabel('$\lambda$', INT,LX);    
    xl.FontSize = xlFS;