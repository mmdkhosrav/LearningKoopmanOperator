%--------------------------------------------------------------------------
% Learnin Koopman Operator
% Ref: Representer Theorem for Learning Koopman Operators 
% Link: https://arxiv.org/abs/2208.01681
%
% This code is used for CV part of Learn_Koopman_Operator.m which solves the 
% following learning problem 
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
function err = CV_LKO(theta,ZGY_CV,learning_type)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% initializations: 
if strcmp(learning_type, 'fro') || ...
          strcmp(learning_type, 'opr') || ...
            strcmp(learning_type, 'nuc')  
    loglambda_est = theta.loglambda;
else
    rnk_est = theta.rnk;    
end
%--------------------------------------------------------------------------
if strcmp(learning_type, 'fro') || strcmp(learning_type, 'rnk')
    Zt  =  ZGY_CV.Zt;
    KGt =  ZGY_CV.KGt;
    Yt  =  ZGY_CV.Yt;
    
    PGinv_GVt = ZGY_CV.PGinv_GVt;
    KXvPt  = ZGY_CV.KXvPt;
    KPtPv  =  ZGY_CV.KPtPv;
    KXpvPv =  ZGY_CV.KXpvPv;
    
    if strcmp(learning_type, 'fro') 
        lambda_est = exp(loglambda_est);
        At = opt_sol_fro(Zt,KGt,Yt,lambda_est);
    else
        [At,~,~] = opt_sol_rnk(Zt,KGt,Yt,rnk_est);
    end
    
    Ct = PGinv_GVt * At;
    E_fro = KXvPt * Ct * KPtPv - KXpvPv;
    err = sum(E_fro(:).^2);

elseif strcmp(learning_type, 'nuc') || strcmp(learning_type, 'opr')
    Ztr  =  ZGY_CV.Ztr;
    KGtr =  ZGY_CV.KGtr;
    Yt   =  ZGY_CV.Yt;
    
    PGinv_GVt = ZGY_CV.PGinv_GVt;
    KXvPt  = ZGY_CV.KXvPt;
    KPtPv  =  ZGY_CV.KPtPv;
    KXpvPv =  ZGY_CV.KXpvPv;
    
    lambda_est = exp(loglambda_est);
    B0 = (eye(size(Ztr))/Ztr) * Yt * (eye(size(KGtr))/KGtr);
    
    if strcmp(learning_type, 'nuc') 
        [At,~,~,~,~,~,~] = opt_sol_nuc(Ztr,KGtr,Yt,lambda_est,B0);
    else
        [At,~,~,~,~,~,~] = opt_sol_opr(Ztr,KGtr,Yt,lambda_est,B0);
    end
    Ct = PGinv_GVt * At;
    E_nuc = KXvPt * Ct * KPtPv - KXpvPv;
    err = sum(E_nuc(:).^2);
end

end
