%--------------------------------------------------------------------------
% Learnin Koopman Operator
% Ref: Representer Theorem for Learning Koopman Operators 
% Link: https://arxiv.org/abs/2208.01681
%
% Example: performance comparison for learning Koopman operator for the EDMD 
% method and the case with Frobenius regularization (Example 4 in ref).
%
% In this example, we consider the PDE describing u(xi,t) as follows
%       du/dt = du/dxi + 0.1 * d^2u/dxi^2,  
% and formulate the learning problem as 
%       min_K E(K) + lambda * R(K),
% where E(.) is the sum squared error loss function and R(.) is the square 
% of Frobenius norm of K, i.e.,  R(K) = ||K||_F^2. We consider two cases, 
% one is lambda = 0, which corresponds to the EDMD method and one with tunned 
% lambda using CV.
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
disp('Example: performance comparison for learning Koopman operator for  ') 
disp('the EDMD method and the case with Frobenius regularization (Example') 
disp('4 in the reference above).                                                         ')
disp(' ')
disp('In this example, we consider the PDE describing u(xi,t) as follows ')
disp('       du/dt = du/dxi + 0.1 * d^2u/dxi^2,                          ')
disp('and formulate the learning problem as                              ')
disp('       min_K      E(K) + lambda * R(K),                            ')
disp('where E(.) is the sum squared error loss function, and R(.) is the ')
disp('square of Frobenius norm of K, i.e.,  R(K) = ||K||_F^2. We consider') 
disp('two cases,  one is lambda = 0, which corresponds to the EDMD method') 
disp('and one with tunned lambda using CV.                               ')
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
%--------------------------------------------------------------------------
% PDE:       du/dt = du/dxi + 0.1 * d^2u/dxi^2  

% discretizing PDE
Dxi = 1e-2;
Dt = 1e-4;

% Nx is the dimension discretized variable x in x_{k+1} = F(x_{k})
Dxi = 1/floor(1/Dxi);  
nx = floor(1/Dxi)+1; 

% discrete first derivative in xi-domain
Dmx = -[eye(nx) zeros(nx,1)] + [zeros(nx,1) eye(nx)];
Dmx = Dmx(:,1:end-1);

% discrete second derivative in xi-domain
D2mx = [eye(nx) zeros(nx,2)] + [zeros(nx,2) eye(nx)];
D2mx = D2mx(:,2:end-1)-2*eye(nx);
D2mx(1,1) = -1; 
D2mx(nx,nx) = -1;

% the matrix for the discretized dynamics
Amx = eye(nx) + Dt/Dxi * Dmx +  0.1 * Dt/Dxi^2 * D2mx;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Generating tajectory of the system for training and validation

% time range
T  = 0: Dt: 1.25;

% initial conditions:
Xi0 = sin(pi*(0:Dxi:1));

% trajectory data
X  = [];   % X:   [x(0) x(1) x(2) ... x(nS-1)      ] (for single trajectory)
Xp = [];   % Xp:  [     x(1) x(2) ... x(nS-1) x(nS)] (for single trajectory)
XXp = [];  % XXp: [x(0) x(1) x(2) ... x(nS-1) x(nS)] (for single trajectory)

% some index set for the trjaectories
idx_X = [];     idx_Xp = [];    idx_X_01 = [];  idx_Xp_01 = []; idx_XXp = [];

for i = 1: size(Xi0,1) % number of trajectories
    
    % initiating dynamical system
    Xi = Xi0(i,:);
    
    idx_X_j = [];       
    idx_Xp_j = [];
    
    % generating the system i-th trajectory
    for j=1:length(T)-1
        % next sample
        Xi = [Xi; DiffE_PDE(Xi(end,:),Amx)];

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

% We have:
% XXp(idx_X_loc,:)  = X  or X  = XXp(idx_X_loc,:) 
% XXp(idx_Xp_loc,:) = Xp or Xp = XXp(idx_Xp_loc,:) 
idx_X_loc  = find(idx_X_01  == 1);
idx_Xp_loc = find(idx_Xp_01 == 1);

% to save data points:
% save('data_DiffE_PDE.mat',...
%     'Xi0','nt', ...
%     'X','Xp','XXp',...
%     'idx_X','idx_Xp',...
%     'idx_X_01','idx_Xp_01','idx_XXp',...
%     'idx_X_loc','idx_Xp_loc')

%--------------------------------------------------------------------------
% % Plotting initial solution
% [T_range,XI_range] = meshgrid(T,0:Dxi:1);
% figure()
% mesh(T_range,XI_range,XXp')
% getframe;

%--------------------------------------------------------------------------
disp('Generating the training and validation trajectory is finished!')
format shortg
c = clock;
disp(['End time: ',num2str(c(1)),'-',num2str(c(2))','-',num2str(c(3)),' ',...
    num2str(c(4)),':',num2str(c(5))])
disp(' ')

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Generating tajectory of the system for test

% initial conditions:
Wi0 = 1 - exp(-(0:Dxi:1));

% trajectory data (similar to above but for test)
W  = [];	Wp = [];	WWp = [];

% some index set for the trjaectories (similar to above but for test)
idx_W = [];     idx_Wp = [];	idx_WWp = [];	idx_W_01 = [];	idx_Wp_01 = [];

for i = 1:size(Wi0,1)
    
    % initiating dynamical system
    Wi = Wi0(i,:);
    
    idx_W_j = [];       
    idx_Wp_j = [];
    
    % generating the system i-th trajectory
    for j=1:length(T)-1
        % next sample
        Wi = [Wi; DiffE_PDE(Wi(end,:),Amx)];

        % collecting the data point        
        idx_W_j  = [idx_W_j  ;  i    j-1];
        idx_Wp_j = [idx_Wp_j ;  i    j];
    end
    
    % concatenate data points of trajectories
    W   = [W;  Wi(1:end-1,:)];
    Wp  = [Wp; Wi(2:end,:)];
    WWp = [WWp; Wi(1:end,:)];
    
    % updating index sets
    idx_W   = [idx_W;    idx_W_j];
    idx_Wp  = [idx_Wp;   idx_Wp_j];
    idx_WWp = [idx_WWp;  idx_W_j;   i   j];    
    idx_W_01   = [idx_W_01;   ones(j,1);    0           ];
    idx_Wp_01  = [idx_Wp_01;  0;            ones(j,1)   ];
end

% We have:
% WWp(idx_W_loc,:)  = W  or W  = WWp(idx_W_loc,:) 
% WWp(idx_Wp_loc,:) = Wp or Wp = WWp(idx_Wp_loc,:)
idx_W_loc  = find(idx_W_01  == 1);
idx_Wp_loc = find(idx_Wp_01 == 1);

%--------------------------------------------------------------------------
% % Plotting test solution
% [T_range,W_range] = meshgrid(T,0:Dxi:1);
% figure()
% mesh(T_range,W_range,WWp')
% getframe;

%--------------------------------------------------------------------------
disp('Generating the test trajectory is finished!')
format shortg
c = clock;
disp(['End time: ',num2str(c(1)),'-',num2str(c(2))','-',num2str(c(3)),' ',...
    num2str(c(4)),':',num2str(c(5))])
disp(' ')

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Choosing kernel type and parameters
kernel_type = 'L';
theta_k = [1 zeros(1,nx)];     

%--------------------------------------------------------------------------
% The observales are defined as 
%       g_l(.) = k(p_l, .),     l = 1,..., nG,
% where k is the kernel function. The matrix P collects p_1,...,p_nG, i.e.,
% P = [p_1',p_2',...,p_nG']'. Here, we set nG = 2 * nX and draw them
% randomly from standard normal distribution in R^nx, i.e., N(0_nx,I_nx).
P = mvnrnd(zeros(1,nx),eye(nx),2*nx);   % data002.mat

% KG is the Gram matrix of p_1,...,p_nG, i.e., KG_(i,j) = k(p_i,p_j)    
KG  = kernel_fun(P, P, theta_k, kernel_type);

% making sure that G is positive definite and generating its square root Gr.
eps_KG = 1e-8; 
[KG,KGr] = mxPD(KG,eps_KG);
nG = size(KG,1);

% choosing observables for training and validation
CVg_ratio = 0.667;  % 2/3 for training, 1/3 for validation 
CVg_step = floor(1/(1 - CVg_ratio));

% some index set for the training and validation observables
idx_P = (1:size(KG,1))';
idx_Pt_loc = idx_P(mod(idx_P,CVg_step) ~= 0);  % training
idx_Pv_loc = idx_P(mod(idx_P,CVg_step) == 0);  % validation
Pt = P(idx_Pt_loc,:);       Pv = P(idx_Pv_loc,:);

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

% Base DATA ---------------------------------------------------------------
% Calculating matrix KXXpP with entries KXXpP_ij = k(x(i),p_j), forci = 0, 
% ... , nS and  j = 1, ... , nG
KXXpP = kernel_fun(XXp, P, theta_k, kernel_type); 

% According to definition of Y, we can take it as a sub-matrix of KXXpP
Y  = KXXpP(idx_Xp_loc,:);

% % trajectory data
% X  = KXXp(idx_X_loc,:);
% Xp = KXXp(idx_Xp_loc,:);

% Calculating matrix KPG with entries KPG_ij = k(x(i),p_j), for i = 0, ... 
% , nS and  j = 1, ... , nG. Note that KPG is a sub-matrix of KXXpP
KPG = kernel_fun(X,P,theta_k,kernel_type);
% KPG = KXXpP(idx_X_loc,:);  % this is an alternative faster way

% Here, we calculate matrices PGinv_GV = inv(G)*KPG and GV = KPG*inv(G)*KPG. 
% For more details, please see Remark 4 in the reference.
PGinv_GV = (eye(size(KG))/KG) * (KPG');
GV = KPG * PGinv_GV;

% To unify the notations, we use Z.
Z = GV; 

% making sure that Z is positive definite and generating its square matrix 
% root Zr. This also improves numerical performance of the implementations.
eps_Z = 1e-8;
[Z,Zr] = mxPD(Z,eps_Z);
nZ = size(Z,1);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% splitting trajectory data points for training and validation
CVx_ratio = 0.85; % 85 percent for training, 15 percent for validation 
CVx_step = floor(1/(1 - CVx_ratio));
 
% some index set for the training and validation trajectory data points
idx_X_vec = (1:size(X,1))';
idx_XXpt_loc = idx_X_vec(mod(idx_X_vec,CVx_step) ~= 0);  % training
idx_XXpv_loc = idx_X_vec(mod(idx_X_vec,CVx_step) == 0);  % validation
Xt  = X(idx_XXpt_loc,:);
Xpt = Xp(idx_XXpt_loc,:);
Xv  = X(idx_XXpv_loc,:); 
Xpv = Xp(idx_XXpv_loc,:);
nst = size(Xt,1);
nsv = size(Xv,1);

% training part of Y
% Yt = kernel_fun(Xpt,Pt,theta_k,kernel_type);
Yt = Y(idx_XXpt_loc,idx_Pt_loc);

% KPGt = kernel_fun(Xt,Pt,theta_k,kernel_type);
KPGt = KPG(idx_XXpt_loc,idx_Pt_loc);

% Calculation of GV for training data
PGinv_GVt = (eye(size(KGt))/KGt) * (KPGt');
GVt = KPGt * PGinv_GVt;

% Calculation of Z matrix for training data
Zt = GVt; 

% making sure that Zt is positive definite and generating its square root 
% Ztr.
eps_Zt = 1e-8;
[Zt,Ztr] = mxPD(Zt,eps_Zt);
nZt = size(Zt,1);

% This matrix is required for calculating validation error 
KXvPt  = kernel_fun(Xv,Pt,theta_k,kernel_type);
% XvPt = PG(idx_XXpv_loc,idx_Pt_loc); % this is an alternative faster way

% This matrix is required for calculating validation error
KXpvPv = kernel_fun(Xpv,Pv,theta_k,kernel_type);
% % KXpP = kernel_fun(Xp,P,theta_k,kernel_type);
% % KXpP = KXXpP(idx_Xp_loc,:);  % this is an alternative faster way
% % KXpvPv = KXpP(idx_XXpv_loc,idx_Pv_loc); % this is an alternative faster way

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% storing all matrices in ZGY
ZGY.X   =  X;               ZGY.Xt   =  Xt;
ZGY.Y   =  Y;               ZGY.Yt   =  Yt;
ZGY.Z   =  Z;               ZGY.Zt   =  Zt;
ZGY.Zr  =  Zr;              ZGY.Ztr  =  Ztr;
ZGY.KG  =  KG;              ZGY.KGt  =  KGt;
ZGY.KGr =  KGr;             ZGY.KGtr =  KGtr;
ZGY.PGinv_GV  =  PGinv_GV;  ZGY.PGinv_GVt = PGinv_GVt;
ZGY.KXvPt  = KXvPt;         ZGY.KPtPv  =  KPtPv;
ZGY.KXpvPv = KXpvPv;

disp('-------------------------------------------------------------------')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Extended dynamic mode decomposition (EDMD):
% We would like to solve the following CONVEX programm
%       min_B ||Z * B * G - Y||_F^2,
% which is the case of Frobenius norm with lambda = 0.

opt.learning_type = 'edmd';
disp('Running EDMD method ...')
format shortg
c = clock;
disp(['Start time: ',num2str(c(1)),'-',num2str(c(2))','-',num2str(c(3)),' ',...
    num2str(c(4)),':',num2str(c(5))])
disp(' ')

[A_edmd, C_edmd] = Learn_Koopman_Operator(ZGY,opt);


disp(' ')
disp('Running EDMD method is finished!')
format shortg
c = clock;
disp(['End time: ',num2str(c(1)),'-',num2str(c(2))','-',num2str(c(3)),' ',...
    num2str(c(4)),':',num2str(c(5))])
disp(' ')
disp('-------------------------------------------------------------------')

%--------------------------------------------------------------------------
% Frobenius norm regularization:
% We would like to solve the following CONVEX programm
%       min_B ||Z * B * G - Y||_F^2 + lambda * ||Zr * B * Gr||_F^2,
% which has a closed form solution. This implementation is highly scalable
% in terms of time and memory. Also, it has a great performance.

opt.learning_type = 'fro';
opt.loglambda_range = [-20 22];
opt.nB0 = 30;

disp('Running with Frobenius norm regularization ...')
disp('This will take a while (several hours on a standard PC). Please be ')
disp('patient.                                                           ') 

format shortg
c = clock;
disp(['Start time: ',num2str(c(1)),'-',num2str(c(2))','-',num2str(c(3)),' ',...
    num2str(c(4)),':',num2str(c(5))])
disp(' ')
[A_fro,C_fro]  = Learn_Koopman_Operator(ZGY,opt);

disp(' ')
disp('Running with Frobenius norm regularization is finished!')
format shortg
c = clock;
disp(['End time: ',num2str(c(1)),'-',num2str(c(2))','-',num2str(c(3)),' ',...
    num2str(c(4)),':',num2str(c(5))])
disp(' ')
disp('-------------------------------------------------------------------')

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% reconstructing W using EDMD estimation of Koopman operator
W_hat_edmd = zeros(size(WWp)); 

W_hat_edmd(1,:) = Wi0; 
gx0CGk = (theta_k(1) + P*(Wi0'))';

for k = 2: length(T)
    gx0CGk = gx0CGk * C_edmd * KG;
    tmp = P\(gx0CGk'-1);
    W_hat_edmd(k,:) = tmp';     
end

Err = W_hat_edmd'-WWp';
L2err_edmd = sum(Err(:).^2)*Dt*Dxi;
disp(['The L2 error for the EDMD method is ', num2str(L2err_edmd)])

%--------------------------------------------------------------------------
% reconstructing W using estimation of Koopman operator with Frobenius norm
% regularization

W_hat_fro = zeros(size(WWp)); 

W_hat_fro(1,:) = Wi0; 
gx0CGk = (theta_k(1) + P*(Wi0'))';

for k = 2: length(T)
    gx0CGk = gx0CGk * C_fro * KG;
    tmp = P\(gx0CGk'-1);
    W_hat_fro(k,:) = tmp';     
end

Err = W_hat_fro'-WWp';
L2err_fro = sum(Err(:).^2) * Dt * Dxi;
disp(['The L2 is for regularization with squared of Frobenius norm is ', ...
    num2str(L2err_fro)])

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% plotting the results

% Font sizes:
tFS   = 15;    % title font sizes
xlFS  = 12;    % xlable font sizes
ylFS  = 12;    % ylable font sizes
xtlFS = 11;    % x-tick font sizes
ytlFS = 11;    % y-tick font sizes

% the range of results data
% max_data = max(max([WWp W_hat_fro W_hat_edmd])); % W_hat_edmd can be +Inf
% min_data = min(min([WWp W_hat_fro W_hat_edmd])); % W_hat_edmd can be -Inf
max_data = max(max([WWp W_hat_fro]));
min_data = min(min([WWp W_hat_fro]));

[T_range,W_range] = meshgrid(T,0:Dxi:1);

N_plot = 33;
f = figure();
f.Position = [250 250 900 350]; 
    subplot(1,N_plot,1:9); 
    surf(T_range,W_range,WWp')
    axis tight; box on; shading interp;
    
    caxis manual
    caxis([min_data max_data]);    
    view(2);
    xlim([0 1.25]);
    ylim([0  1]);
    ax = gca;
    ax.TickLabelInterpreter = LX;
    ax.XTick = [0 0.5 1];
    ax.XTickLabel(1) = {'$0$'};    
    ax.XTickLabel(2) = {'$0.5$'};    
    ax.XTickLabel(3) = {'$1$'}; 
    xtl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',xtl,'fontsize',xtlFS,'FontName','Times New Roman');
    ax.YTick = [0 0.5 1];
    ax.YTickLabel(1) = {'$0$'};    
    ax.YTickLabel(2) = {'$0.5$'};    
    ax.YTickLabel(3) = {'$1$'}; 
    ytl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',ytl,'fontsize',ytlFS,'FontName','Times New Roman');
    xl = xlabel('$t$', INT,LX);
    xl.FontSize = xlFS;  
    yl = ylabel('$\xi$', INT,LX);
    yl.FontSize = ylFS;   
    t = title('$u(\xi,t)$', INT,LX);
    t.FontSize = tFS;
    
    %----------------------------------------------------------------------
    subplot(1,N_plot,11:20); 
    surf(T_range,W_range,W_hat_edmd')
    axis tight; box on; shading interp;
    
    caxis manual
    caxis([min_data max_data]);
    view(2);
    xlim([0 1.25]);
    ylim([0  1]);
    ax = gca;
    ax.TickLabelInterpreter = LX;
    ax.XTick = [0 0.5 1];
    ax.XTickLabel(1) = {'$0$'};    
    ax.XTickLabel(2) = {'$0.5$'};    
    ax.XTickLabel(3) = {'$1$'}; 
    xtl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',xtl,'fontsize',xtlFS,'FontName','Times New Roman');
    ax.YTick = [0 0.5 1];
    ax.YTickLabel(1) = {''};    
    ax.YTickLabel(2) = {''};    
    ax.YTickLabel(3) = {''}; 
    ytl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',ytl,'fontsize',ytlFS,'FontName','Times New Roman');
    
    xl = xlabel('$t$', INT,LX);
    xl.FontSize = xlFS;  
    % yl = ylabel('$\xi$', INT,LX);
    % yl.FontSize = ylFS;   
    t = title('$\tilde{u}(\xi,t)$', INT,LX);
    t.FontSize = tFS;

    %----------------------------------------------------------------------
    subplot(1,N_plot,22:N_plot);  
    surf(T_range,W_range,W_hat_fro');
    axis tight; box on; shading interp;
    
    caxis manual
    caxis([min_data max_data]);
    colorbar;
    view(2);
    xlim([0 1.25]);
    ylim([0  1]);
    ax = gca;
    ax.TickLabelInterpreter = LX;
    ax.XTick = [0 0.5 1];
    ax.XTickLabel(1) = {'$0$'};    
    ax.XTickLabel(2) = {'$0.5$'};    
    ax.XTickLabel(3) = {'$1$'}; 
    xtl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',xtl,'fontsize',xtlFS,'FontName','Times New Roman');
    ax.YTick = [0 0.5 1];
    ax.YTickLabel(1) = {''};    
    ax.YTickLabel(2) = {''};    
    ax.YTickLabel(3) = {''}; 
    ytl = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',ytl,'fontsize',ytlFS,'FontName','Times New Roman');
    xl = xlabel('$t$', INT,LX);
    xl.FontSize = xlFS;      
    % yl = ylabel('$\xi$', INT,LX);
    % yl.FontSize = ylFS;       
    t = title('$\hat{u}(\xi,t)$', INT,LX);
    t.FontSize = tFS;

