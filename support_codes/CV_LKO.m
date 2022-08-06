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
function err = CV_LKO(theta,ZGY,learning_type)
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
    if strcmp(learning_type, 'fro') 
        lambda_est = exp(loglambda_est);
        At = opt_sol_fro(ZGY.Zt,ZGY.KGt,ZGY.Yt,lambda_est);
    else
        [At,~,~] = opt_sol_rnk(ZGY.Zt,ZGY.KGt,ZGY.Yt,rnk_est);
    end
    
    Ct = (ZGY.PGinv_GVt) * At;
    E_val = (ZGY.KXvPt) * Ct * (ZGY.KPtPv) - (ZGY.KXpvPv);
    err = sum(E_val(:).^2);

elseif strcmp(learning_type, 'nuc') || strcmp(learning_type, 'opr')
    lambda_est = exp(loglambda_est);
    B0 = (eye(size(ZGY.Ztr))/(ZGY.Ztr)) * (ZGY.Yt) * ...
                (eye(size(ZGY.KGtr))/(ZGY.KGtr));
    
    if strcmp(learning_type, 'nuc') 
        [At,~,~,~,~,~,~] = opt_sol_nuc(ZGY.Ztr,ZGY.KGtr,ZGY.Yt,lambda_est,B0);
    else
        [At,~,~,~,~,~,~] = opt_sol_opr(ZGY.Ztr,ZGY.KGtr,ZGY.Yt,lambda_est,B0);
    end
    Ct = (ZGY.PGinv_GVt) * At;
    E_val = (ZGY.KXvPt) * Ct * (ZGY.KPtPv) - (ZGY.KXpvPv);
    err = sum(E_val(:).^2);
end

end
