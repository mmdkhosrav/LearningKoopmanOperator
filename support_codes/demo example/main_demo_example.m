%--------------------------------------------------------------------------
% Learnin Koopman Operator
% Ref: Representer Theorem for Learning Koopman Operators 
% Link: https://arxiv.org/abs/2208.01681
%
% This code is a demo example for solving the following learning problem 
%   min_K   E(K) + lambda * R(K),
%   s.t.    K in C,
% where E(.) is the sum squared error loss function, C is a constraint R(.)
% the regularization term. For C and R, we have the following cases:
%   1. R: square of operator norm of K, i.e.,  R(K) = ||K||^2 
%   2. R: square of Frobenius norm of K, i.e.,  R(K) = ||K||_F^2 
%   3. R: nuclear norm of K, i.e.,  R(K) = ||K||_* 
%   4. C: rank of K, i.e.,  C = {K | rank(K) <= r} 
%   5. R&C: if R=0 and C = L(K) (the whole space of bounded operators), the
%      the program is equivalent to EDMD approach
%
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------

close all;  clear all;  clc;

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
addpath(currentFolder(1:end-12))

%--------------------------------------------------------------------------
% SNR level for adding noise to the system trajectories:  
dBvalx = 10;
rng(1000 * dBvalx); % for repeatability

%--------------------------------------------------------------------------
% Generating tajectory of the system
% range:
x1_min = -3;        x1_max = 3;
x2_min = -2.5;      x2_max = 2.5;

eps_X_lim = 1e-2; 

% initial point:
Xi0 = [ -2     2; 
         0    -1];
     
% number of trajectories     
nt = size(Xi0,1);

% sampling time for discretizing trajectory
dT = 0.2;

% time range
T  = 0:dT:10;


% trajectory data
X  = [];   % X:   [x(0) x(1) x(2) ... x(nS-1)      ] (for single trajectory)
Xp = [];   % Xp:  [     x(1) x(2) ... x(nS-1) x(nS)] (for single trajectory)
XXp = [];  % XXp: [x(0) x(1) x(2) ... x(nS-1) x(nS)] (for single trajectory)

% some index set for the trjaectories
idx_X = [];     idx_Xp = [];    idx_X_01 = [];  idx_Xp_01 = []; idx_XXp = [];


for i = 1:nt
    % initiating dynamical system
    Xi = Xi0(i,:);
    
    % some index set
    idx_X_j = [];       
    idx_Xp_j = [];

    % generating the system i-th trajectory
    for j=1:length(T)-1
        % next sample
        Xi = [Xi; DiffE_VP(Xi(end,:),dT)];

        % In case we want to avoid similar trajectory point, we can set         
        % eps_X to small value and use following lines:
        % eps_X = norm(Xi(end-1,:)-Xi(end,:));
        % if eps_X < eps_X_lim
        %     break;
        % end

        % collecting the data point
        idx_X_j  = [idx_X_j  ;  i    j-1];
        idx_Xp_j = [idx_Xp_j ;  i    j];
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

idx_X_loc  = find(idx_X_01  == 1);
idx_Xp_loc = find(idx_Xp_01 == 1);

% We have:
% XXp(idx_X_loc,:)  = X  or X  = XXp(idx_X_loc,:) 
% XXp(idx_Xp_loc,:) = Xp or Xp = XXp(idx_Xp_loc,:) 

%--------------------------------------------------------------------------
% Choosing kernel type and parameters
kernel_type = 'M5';    % MAtern kernel 5/2
theta_k = [1 1 0 0 0];

%--------------------------------------------------------------------------
% The observales are defined as 
%       g_l(.) = k(p_l, .),     l = 1,..., nG,
% where k is the kernel function. The matrix P collects p_1,...,p_nG, i.e.,
% P = [p_1',p_2',...,p_nG']'. We pick P as a grid in a suitable range.

% discretizing the x1-axis and x2-axis
Nx1_P = 10;     Nx2_P = 10;   
% deciding on the range
x1_P_min = -2.5;              x1_P_max = 2.5;
x2_P_min = -2.5;              x2_P_max = 2.5;
tmp_x1_P = (0:1/Nx1_P:1)';          tmp_x2_P = (0:1/Nx2_P:1)'; 
x1_P_range = ( - x1_P_min + x1_P_max) * tmp_x1_P + x1_P_min;
x2_P_range = ( - x2_P_min + x2_P_max) * tmp_x2_P + x2_P_min;
[X1_P, X2_P] = meshgrid(x1_P_range,x2_P_range);
P = [X1_P(:) X2_P(:)];

% KG is the Gram matrix of p_1,...,p_nG, i.e., KG_(i,j) = k(p_i,p_j)  
KG  = kernel_fun(  P, P, theta_k, kernel_type);

% making sure that KG is positive definite and generating its square root 
% matrix KGr.
eps_KG = 1e-8; 
[KG,KGr] = mxPD(KG,eps_KG);
nG = size(KG,1);

% choosing observables for training and validation
CVg_ratio = 0.50;      % 50 percent for training, 50 percent for validation 
CVg_step = floor(1/(1 - CVg_ratio) - 1e-8)+ 1;

% some index set for the training and validation observables
idx_P = (1:size(KG,1))';
idx_Pt_loc = idx_P(mod(idx_P,CVg_step) ~= 0);  % training
idx_Pv_loc = idx_P(mod(idx_P,CVg_step) == 0);  % validation
Pt = P(idx_Pt_loc,:);   % observables for training   
Pv = P(idx_Pv_loc,:);   % observables for validation   

KGt  = kernel_fun(Pt,Pt,theta_k,kernel_type);
% KGt  = KG(idx_Pt_loc,idx_Pt_loc);
eps_KGt = 1e-8;
[KGt,KGtr] = mxPD(KGt,eps_KGt);
nGt = size(KGt,1);

KPtPv  = kernel_fun(Pt,Pv,theta_k,kernel_type);
% KPtPv = KG(idx_Pt_loc,idx_Pv_loc);

%--------------------------------------------------------------------------
% observables for test
Nx1_Ptst = 50;      Nx2_Ptst = 50;      s_Ptst = 1;

x1_Ptst_min = x1_min * s_Ptst;      x1_Ptst_max = x1_max * s_Ptst;
x2_Ptst_min = x2_min * s_Ptst;      x2_Ptst_max = x2_max * s_Ptst;
tmp_x1_Ptst = (0:1/Nx1_Ptst:1)';    tmp_x2_Ptst = (0:1/Nx2_Ptst:1)'; 
x1_Ptst_range = ( - x1_Ptst_min + x1_Ptst_max) * tmp_x1_Ptst + x1_Ptst_min;
x2_Ptst_range = ( - x2_Ptst_min + x2_Ptst_max) * tmp_x2_Ptst + x2_Ptst_min;

[X1_Ptst, X2_Ptst] = meshgrid(x1_Ptst_range,x2_Ptst_range);

Ptst = [X1_Ptst(:) X2_Ptst(:)];

KGtst  = kernel_fun(  Ptst, Ptst, theta_k, kernel_type);
eps_KGtst = 1e-8; 
[KGtst,KGtstr] = mxPD(KGtst,eps_KGtst);

%--------------------------------------------------------------------------
% XGrid for test: Xtst
Nx1_Xtst = 100;    Nx2_Xtst = 100;    s_Xtst = 1;

x1_Xtst_min = x1_min * s_Xtst;     x1_Xtst_max = x1_max * s_Xtst;
x2_Xtst_min = x2_min * s_Xtst;     x2_Xtst_max = x2_max * s_Xtst;
tmp_x1_Xtst = (0:1/Nx1_Xtst:1)';   tmp_x2_Xtst = (0:1/Nx2_Xtst:1)'; 
dx1dx2_Xtst = (-x1_Xtst_min + x1_Xtst_max)/Nx1_Xtst ...
                * (-x2_Xtst_min + x2_Xtst_max)/Nx2_Xtst;   

x1_Xtst_range = ( - x1_Xtst_min + x1_Xtst_max) * tmp_x1_Xtst + x1_Xtst_min;
x2_Xtst_range = ( - x2_Xtst_min + x2_Xtst_max) * tmp_x2_Xtst + x2_Xtst_min;

[X1_Xtst, X2_Xtst] = meshgrid(x1_Xtst_range,x2_Xtst_range);

Xtst  = [X1_Xtst(:) X2_Xtst(:)];

Xptst = DiffE_VP(Xtst,dT);

KXptstPtst  = kernel_fun(Xptst, Ptst, theta_k, kernel_type);
KXtstP      = kernel_fun(Xtst, P, theta_k, kernel_type);
KPPtst      = kernel_fun(P, Ptst, theta_k, kernel_type); 
%==========================================================================
%--------------------------------------------------------------------------
% adding noise to XXp
XXp_NL = XXp; 
W_mx = randn(size(XXp));
W_mx = 10^((snr(XXp(:),W_mx(:))-dBvalx)/20) * W_mx;
W_mx_var = std(W_mx(:))^2;
XXp = XXp + W_mx;                       % noisy output
SNR_XXp = snr(XXp(:),W_mx(:));

% Base DATA ---------------------------------------------------------------

XXpP = kernel_fun(XXp, P, theta_k, kernel_type); %XY

X  = XXp(idx_X_loc,:);
Xp = XXp(idx_Xp_loc,:);
Y  = XXpP(idx_Xp_loc,:);

PG = kernel_fun(X,P,theta_k,kernel_type);
% PG = XXpP(idx_X_loc,:);

XpP = kernel_fun(Xp,P,theta_k,kernel_type);
% XpP = XXpP(idx_Xp_loc,:);

% We know that GV = PG * inv(KG) * PG';
PGinv_GV = (eye(size(KG))/KG) * (PG');
GV = PG * PGinv_GV;

Z = GV; 
eps_Z = 1e-8;
[Z,Zr] = mxPD(Z,eps_Z);
nZ = size(Z,1);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% CV DATA

CVx_ratio = 0.75;
CVx_step = floor(1/(1 - CVx_ratio));
 
idx_X_vec = (1:size(X,1))';
idx_XXpt_loc = idx_X_vec(mod(idx_X_vec,CVx_step) ~= 0);  % training
idx_XXpv_loc = idx_X_vec(mod(idx_X_vec,CVx_step) == 0);  % validation
Xt  = X(idx_XXpt_loc,:);
Xpt = Xp(idx_XXpt_loc,:);
Xv  = X(idx_XXpv_loc,:); 
Xpv = Xp(idx_XXpv_loc,:);

nst = size(Xt,1);
nsv = size(Xv,1);

%--------------------------------------------------------------------------
% Yt = kernel_fun(Xpt,Pt,theta_k,kernel_type);
Yt = Y(idx_XXpt_loc,idx_Pt_loc);

%--------------------------------------------------------------------------
% PGt = kernel_fun(Xt,Pt,theta_k,kernel_type);
PGt = PG(idx_XXpt_loc,idx_Pt_loc);

PGinv_GVt = (eye(size(KGt))/KGt) * (PGt');
GVt = PGt * PGinv_GVt;

Zt = GVt; 
eps_Zt = 1e-8;
[Zt,Ztr] = mxPD(Zt,eps_Zt);
nZt = size(Zt,1);

% for validation 
KXvPt  = kernel_fun(Xv,Pt,theta_k,kernel_type);
% XvPt = PG(idx_XXpv_loc,idx_Pt_loc);

KXpvPv = kernel_fun(Xpv,Pv,theta_k,kernel_type);
% XpvPv = XpP(idx_XXpv_loc,idx_Pv_loc);

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
ZGY.X   =  X;               ZGY.Xt   =  Xt;
ZGY.Y   =  Y;               ZGY.Yt   =  Yt;
ZGY.Z   =  Z;               ZGY.Zt   =  Zt;
ZGY.Zr  =  Zr;              ZGY.Ztr  =  Ztr;
ZGY.KG  =  KG;              ZGY.KGt  =  KGt;
ZGY.KGr =  KGr;             ZGY.KGtr =  KGtr;
ZGY.PGinv_GV  =  PGinv_GV;  ZGY.PGinv_GVt = PGinv_GVt;
ZGY.KXvPt  = KXvPt;         ZGY.KPtPv  =  KPtPv;
ZGY.KXpvPv = KXpvPv;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% edmd
opt.learning_type = 'edmd';

[A_edmd, C_edmd] = Learn_Koopman_Operator(ZGY,opt);

EE_edmd = KXtstP * C_edmd * KPPtst - KXptstPtst;
err_MSE_edmd  = sum(EE_edmd(:).^2)* dx1dx2_Xtst/nG;
err_R2_edmd   = 100 * (1- ((mean(EE_edmd(:).^2)).^0.5 / ...
            (mean( (KXptstPtst(:)-mean(KXptstPtst(:))).^2)).^0.5));

disp(opt.learning_type)
%--------------------------------------------------------------------------
opt.learning_type = 'fro';
opt.loglambda_range = [-20 22];
opt.nB0 = 30;

[A_fro,C_fro]  = Learn_Koopman_Operator(ZGY,opt);
EE_fro = KXtstP * C_fro * KPPtst - KXptstPtst;
err_MSE_fro  = sum(EE_fro(:).^2)* dx1dx2_Xtst/nG;
err_R2_fro   = 100 * (1- (mean(EE_fro(:).^2)).^0.5/std(KXptstPtst(:)));

disp(opt.learning_type)
%--------------------------------------------------------------------------
opt.learning_type = 'nuc';
opt.loglambda_range = [-20 22];
opt.nB0 = 30;

[A_nuc,C_nuc]  = Learn_Koopman_Operator(ZGY,opt);
EE_nuc = KXtstP * C_nuc * KPPtst - KXptstPtst;
err_MSE_nuc  = sum(EE_nuc(:).^2)* dx1dx2_Xtst/nG;
err_R2_nuc   = 100 * (1- (mean(EE_nuc(:).^2)).^0.5/std(KXptstPtst(:)));

disp(opt.learning_type)
%--------------------------------------------------------------------------
opt.learning_type = 'opr';
opt.loglambda_range = [-20 22];
opt.nB0 = 30;

[A_opr,C_opr]  = Learn_Koopman_Operator(ZGY,opt);
EE_opr = KXtstP * C_opr * KPPtst - KXptstPtst;
err_MSE_opr  = sum(EE_opr(:).^2)* dx1dx2_Xtst/nG;
err_R2_opr   = 100 * (1- (mean(EE_opr(:).^2)).^0.5/std(KXptstPtst(:)));

disp(opt.learning_type)
%--------------------------------------------------------------------------
opt.learning_type = 'rnk';
opt.rnk = [1 size(KGtr,2)];
opt.nB0 = 30;

[A_rnk,C_rnk]  = Learn_Koopman_Operator(ZGY,opt);
EE_rnk = KXtstP * C_rnk * KPPtst - KXptstPtst;
err_MSE_rnk  = sum(EE_rnk(:).^2)* dx1dx2_Xtst/nG;
err_R2_rnk   = 100 * (1- (mean(EE_rnk(:).^2)).^0.5/std(KXptstPtst(:)));

disp(opt.learning_type)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%==========================================================================
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
fn = ['data_',num2str(dBvalx),'.mat'];
save(fn,...
    'A_edmd','C_edmd','err_MSE_edmd','err_R2_edmd',...
    'A_fro', 'C_fro', 'err_MSE_fro', 'err_R2_fro',...
    'A_nuc', 'C_nuc', 'err_MSE_nuc', 'err_R2_nuc',...
    'A_opr', 'C_opr', 'err_MSE_opr', 'err_R2_opr',...
    'A_rnk', 'C_rnk', 'err_MSE_rnk', 'err_R2_rnk',...
    'W_mx')

IterTime = toc
