function IterTime = ...
    LKR_loop_stbSI(MC_num,...
                     XXp_NL,dBvalx,P,KG,KGr,KXtstP,KPPtst,KXptstPtst,dx1dx2_Xtst,...
                     theta_k, kernel_type,...
                     nG,idx_X_loc,idx_Xp_loc)
                 
%--------------------------------------------------------------------------             
tic;   

rng(1000 * dBvalx + MC_num); % for repeatability

% creating the file name
MC_num_str = num2str(MC_num+1000); MC_num_str = MC_num_str(2:end);
fn = ['data_',MC_num_str,'.mat'];

% adding noise to XXp such that we have the SNR of dBvalx
XXp = XXp_NL; 
W_mx = randn(size(XXp));
W_mx = 10^((snr(XXp(:),W_mx(:))-dBvalx)/20) * W_mx;
W_mx_var = std(W_mx(:))^2;
XXp = XXp + W_mx;                       % noisy output
SNR_XXp = snr(XXp(:),W_mx(:));

% Base DATA ---------------------------------------------------------------
% Calculating matrix KXXpP with entries KXXpP_ij = k(x(i),p_j), forci = 0, 
% ... , nS and  j = 1, ... , nG
KXXpP = kernel_fun(XXp, P, theta_k, kernel_type); 

% According to definition of Y, we can take it as a sub-matrix of KXXpP
Y  = KXXpP(idx_Xp_loc,:);

% trajectory data
X  = XXp(idx_X_loc,:);
Xp = XXp(idx_Xp_loc,:);

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
% Extended dynamic mode decomposition (EDMD):
% We would like to solve the following CONVEX programm
%       min_B ||Z * B * G - Y||_F^2,
% which is the case of Frobenius norm with lambda = 0.

A_edmd = opt_sol_fro(Z,KG,Y,0);
C_edmd = PGinv_GV * A_edmd;

EE_edmd = KXtstP * C_edmd * KPPtst - KXptstPtst;
err_MSE_edmd  = sum(EE_edmd(:).^2)* dx1dx2_Xtst/nG;
err_R2_edmd   = 100 * (1- (mean(EE_edmd(:).^2)).^0.5 /std(KXptstPtst(:)));
eig_edmd = eig(C_edmd * KG);

%--------------------------------------------------------------------------
% Frobenius norm regularization:
% We would like to solve the following CONVEX programm
%       min_B ||Z * B * G - Y||_F^2 + lambda * ||Zr * B * Gr||_F^2,
% which has a closed form solution. This implementation is highly scalable
% in terms of time and memory. Also, it has a great performance.
 
lambda_fro = 1e-6;

A_fro = opt_sol_fro(Z,KG,Y,lambda_fro);
C_fro = PGinv_GV * A_fro;

EE_fro = KXtstP * C_fro * KPPtst - KXptstPtst;
err_MSE_fro  = sum(EE_fro(:).^2)* dx1dx2_Xtst/nG;
err_R2_fro   = 100 * (1- (mean(EE_fro(:).^2)).^0.5/std(KXptstPtst(:)));
eig_fro = eig(C_fro * KG);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Frobenius norm regularization with stability-inducing constraints:
% We would like to solve the following CONVEX programm
%       min_B   ||Zr * B * Gr - Y||_F^2 + lambda * ||B||_F^2,
%       s.t.    ||B|| <= 1 - eps
% which can be solved efficienly.

eps_stb = 1e-5;

% deciding whether solve the problem using cvx + mosek or by fmincon
solving_cvx_and_mosek = 0;

% cvx + mosek
if solving_cvx_and_mosek
    disp('Solving using cvx & mosek')
    cvx_solver mosek
    cvx_begin quiet
        variable B(nZ,nG)    
        minimize(sum_square(vec(Zr * B * KGr - Y)) + lambda_fro * sum_square(vec(B)))
        subject to
        norm(B) <= 1 - eps_stb
    cvx_end
    B_cvx = B;

    A_cvx = (eye(nZ)/Zr) * B * (eye(nG)/KGr);

    C_cvx = PGinv_GV * A_cvx;
    EE_cvx = KXtstP * C_cvx * KPPtst - KXptstPtst;
    err_MSE_cvx  = sum(EE_cvx(:).^2)* dx1dx2_Xtst/nG;
    err_R2_cvx   = 100 * (1- (mean(EE_cvx(:).^2)).^0.5/std(KXptstPtst(:)));
    eig_cvx = eig(C_cvx * KG);
    eig_cvx_max = max(abs(eig_cvx));
    
    % saving the results
    save(fn,...
        'A_edmd','C_edmd','err_MSE_edmd','eig_edmd',...
        'A_fro', 'C_fro', 'err_MSE_fro', 'eig_fro',...
        'A_cvx', 'C_cvx', 'err_MSE_cvx', 'eig_cvx',...
        'W_mx')
% fmincon
else
    % disp('Solving using fmincon')
    options = optimoptions(   'fmincon'...
                            , 'Display','off'... 'off'  'iter'
                            , 'Algorithm','interior-point' ...
                            , 'SpecifyObjectiveGradient',true ...
                            , 'SpecifyConstraintGradient',true ...
                            , 'MaxFunctionEvaluations', 1e5...  
                            , 'MaxIterations', 1e4 ...
                            , 'OptimalityTolerance', 1e-10 ...  
                            , 'ConstraintTolerance', 1e-10 ...
                            , 'FunctionTolerance', 1e-10 ...
                            , 'StepTolerance', 1e-10 ...
                            , 'HessianApproximation', {'lbfgs',15} ...
                          );

    B0 = (eye(nZ)/Zr) * Y * (eye(nG)/KGr);
    x0 = B0(:);
    Aineq = [];
    bineq = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];

    [x,fvalB2,exitflag,output] = ...
            fmincon(@(x)obj_fro_ZrGr(x,Y,Zr,KGr,lambda_fro), x0...
                                ,Aineq,bineq,Aeq,beq,lb,ub,...
                                @(x)cns_opr_ZrGr(x,nZ,nG,eps_stb),...
                                options);

    B_star = reshape(x,nZ,nG);
    A_star = (eye(nZ)/Zr) * B_star * (eye(nG)/KGr);

    C_star = PGinv_GV * A_star;
    EE_star = KXtstP * C_star * KPPtst - KXptstPtst;
    err_MSE_star  = sum(EE_star(:).^2)* dx1dx2_Xtst/nG;
    err_R2_star   = 100 * (1- (mean(EE_star(:).^2)).^0.5/std(KXptstPtst(:)));
    eig_star = eig(C_star * KG);
    eig_star_max = max(abs(eig_star));

    % saving the results
    save(fn,...
        'A_edmd','C_edmd','err_MSE_edmd', 'eig_edmd',...
        'A_fro', 'C_fro', 'err_MSE_fro',  'eig_fro',...
        'A_star','C_star','err_MSE_star', 'eig_star',...
        'W_mx')
end
MC_str = num2str(1000+MC_num);
MC_str = MC_str(2:end);
disp(['iter. ',MC_str  ,' done!'])
IterTime = toc;
end