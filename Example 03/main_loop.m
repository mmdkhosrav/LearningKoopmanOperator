%--------------------------------------------------------------------------
% Learnin Koopman Operator
% Ref: Representer Theorem for Learning Koopman Operators 
% Link: https://arxiv.org/abs/2208.01681
%
% This code concerns evaluating the impact of including side-information on
% the learning of Koopman operators. Therefore, we consider the Nicholson 
% and Bailey model for the host and parasitoids dynamics, which is a stable
% system, and subsequently, we compare the performance for the following 
% learning Koopman operator approaches:                     
%    1. EDMD method,                                                 
%    2. Frobenius regularization,                                    
%    3. Frobenius regularization with stability-inducing constraints.  
% More details are provided in Example 4 in the reference above.  
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
disp('This code concerns the impact of including side-information on the ')
disp('learning of Koopman operators. Therefore, we consider the Nicholson')
disp('and Bailey model for the host and parasitoids dynamics, which is a ')
disp('stable system, and subsequently, compare the performance for the   ')
disp('following learning Koopman operator approaches:                    ') 
disp('   1. EDMD method,                                                 ')
disp('   2. Frobenius regularization,                                    ')
disp('   3. Frobenius regularization with stability-inducing constraints.')                    
disp('More details are provided in Example 4 in the reference above.     ')
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
disp(['Start time: ',num2str(c(1)),'-',num2str(c(2))','-',num2str(c(3)),...
    ' ',...
    num2str(c(4)),':',num2str(c(5))])
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
% adding the path of support files
currentFolder = pwd;
addpath([currentFolder(1:end-10),'support_codes'])

CVx_ratio = 0.51;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Generating tajectory of the system
% range:
x1_min = 0;      x1_max = 0.55;
x2_min = 0;      x2_max = 0.3;

% initial point:
Xi0 = [0.5, 0.05];

% number of trajectories     
nt = size(Xi0,1);

% time range
T  = 0:1:100;

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
        Xi = [Xi; DiffE_HP(Xi(end,:))];
  
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
%--------------------------------------------------------------------------
% Choosing kernel type and parameters
kernel_type = 'SE';
theta_k = [0.1 0 0 0 0];

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% The observales are defined as 
%       g_l(.) = k(p_l, .),     l = 1,..., nG,
% where k is the kernel function. The matrix P collects p_1,...,p_nG, i.e.,
% P = [p_1',p_2',...,p_nG']'. We pick P as a grid in a suitable range.

% discretizing the x1-axis and x2-axis for builing P
Nx1_P = 8;     
Nx2_P = 8;   
% deciding on the range
x1_P_min = 0.30;                    x1_P_max = 0.5;
x2_P_min = 0.05;                    x2_P_max = 0.25;
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
CVg_ratio = 0.50;      % 2/3 for training, 1/3 for validation 
CVg_step = floor(1/(1 - CVg_ratio))+1;

% some index set for the training and validation observables
idx_P = (1:size(KG,1))';
idx_Pt_loc = idx_P(mod(idx_P,CVg_step) ~= 0);  % training
idx_Pv_loc = idx_P(mod(idx_P,CVg_step) == 0);  % validation
Pt = P(idx_Pt_loc,:);   % observables for training   
Pv = P(idx_Pv_loc,:);   % observables for validation   

% KGt is the Gram matrix of training observables    
KGt  = kernel_fun(Pt,Pt,theta_k,kernel_type);
% KGt  = KG(idx_Pt_loc,idx_Pt_loc); % this is an alternative faster way

% making sure that KGt is positive definite and generating its square root 
% KGtr. 
eps_KGt = 1e-8;
[KGt,KGtr] = mxPD(KGt,eps_KGt);
nGt = size(KGt,1);

% cross inner products of training and validatin observables  
KPtPv  = kernel_fun(Pt,Pv,theta_k,kernel_type);
% KPtPv = KG(idx_Pt_loc,idx_Pv_loc); % this is an alternative faster way

%--------------------------------------------------------------------------
% We need observables for testing the final outcome. These observables are 
% built similar to the previous ones but with much finer grid. The matrix
% is denoted by Ptst.

% discretizing the x1-axis and x2-axis for builing P
Nx1_Ptst = 20;      
Nx2_Ptst = 20;    
% deciding on the range
x1_Ptst_min = x1_min;               x1_Ptst_max = x1_max;
x2_Ptst_min = x2_min;               x2_Ptst_max = x2_max;
tmp_x1_Ptst = (0:1/Nx1_Ptst:1)';    tmp_x2_Ptst = (0:1/Nx2_Ptst:1)'; 
x1_Ptst_range = ( - x1_Ptst_min + x1_Ptst_max) * tmp_x1_Ptst + x1_Ptst_min;
x2_Ptst_range = ( - x2_Ptst_min + x2_Ptst_max) * tmp_x2_Ptst + x2_Ptst_min;
[X1_Ptst, X2_Ptst] = meshgrid(x1_Ptst_range,x2_Ptst_range);
Ptst = [X1_Ptst(:) X2_Ptst(:)];

% KGtst is the Gram matrix of test observables 
KGtst  = kernel_fun(  Ptst, Ptst, theta_k, kernel_type);

% making sure that KGtst is positive definite and generating its square root 
% KGtstr.
eps_KGtst = 1e-8; 
[KGtst,KGtstr] = mxPD(KGtst,eps_KGtst);
    
%--------------------------------------------------------------------------
% We need an x-domain grid for test through numerical calculation of MSE. 
% This grid is built similar to P and Ptst

% discretizing the x1-axis and x2-axis for builing P
Nx1_Xtst = 50;    
Nx2_Xtst = 50;    
% deciding on the range
x1_Xtst_min = x1_min;              x1_Xtst_max = x1_max;
x2_Xtst_min = x2_min;              x2_Xtst_max = x2_max;
tmp_x1_Xtst = (0:1/Nx1_Xtst:1)';   tmp_x2_Xtst = (0:1/Nx2_Xtst:1)'; 
dx1dx2_Xtst = (-x1_Xtst_min + x1_Xtst_max)/Nx1_Xtst ...
                * (-x2_Xtst_min + x2_Xtst_max)/Nx2_Xtst;   
x1_Xtst_range = ( - x1_Xtst_min + x1_Xtst_max) * tmp_x1_Xtst + x1_Xtst_min;
x2_Xtst_range = ( - x2_Xtst_min + x2_Xtst_max) * tmp_x2_Xtst + x2_Xtst_min;
[X1_Xtst, X2_Xtst] = meshgrid(x1_Xtst_range,x2_Xtst_range);
Xtst  = [X1_Xtst(:) X2_Xtst(:)];

% x-domain grid test shifted by dynamics for one time sample.
Xptst = DiffE_HP(Xtst);

% some kernel inner product calculation to be used later in obtaining MSE
KXptstPtst  = kernel_fun( Xptst , Ptst , theta_k, kernel_type);
KXtstP      = kernel_fun( Xtst  , P    , theta_k, kernel_type);
KPPtst      = kernel_fun( P     , Ptst , theta_k, kernel_type); 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% SNR level for adding noise to the system trajectories:  
dBvalx = 30;
rng(1000 * dBvalx); % for repeatability

% adding noise to XXp such that we have the SNR of dBvalx
XXp_NL = XXp;               % a noiseless backup for XXp    

%--------------------------------------------------------------------------

parfor MC_num = 1:450
    LKR_loop_stbSI(MC_num,...
         XXp_NL,dBvalx,P,KG,KGr,KXtstP,KPPtst,KXptstPtst,dx1dx2_Xtst,...
         theta_k, kernel_type,nG,idx_X_loc,idx_Xp_loc);    
end
