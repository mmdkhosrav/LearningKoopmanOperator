%--------------------------------------------------------------------------
% Learnin Koopman Operator
% Ref: Representer Theorem for Learning Koopman Operators 
% Link: https://arxiv.org/abs/2208.01681
%
% This code solves the following learning problem 
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
% This code will be improved in future!
function  [A,C]  = Learn_Koopman_Operator(ZGY,opt)
%--------------------------------------------------------------------------
% ZGY:  the variable which contains all data and preprocessings, including 
%           Z,G,Y,Zr,Gr,... for mor details see demo example.  
% opt:  the variable which contains options including
%           opt.learning_type:    this can be either 'fro' (this is default
%                      & also recommended), 'nuc', 'opr','rnk' and 'edmd'.                            
%           opt.loglambda_range:  range of log(lambda), where lambda is the
%                      weight of regularization.
%           opt.rnk:   min and max of rank range, for the case with rank 
%                      constraint
%           opt.nBO:    number of Bayesian optimization iteration for CV
%--------------------------------------------------------------------------
% detecting learning type problem in terms of regularization and rank
if ~isfield(opt, 'learning_type')
    learning_type = 'fro';
elseif strcmp(opt.learning_type, 'fro')
    learning_type = 'fro';
elseif strcmp(opt.learning_type, 'nuc') 
    learning_type = 'nuc';
elseif strcmp(opt.learning_type, 'opr') 
    learning_type = 'opr';
elseif strcmp(opt.learning_type, 'rnk') 
    learning_type = 'rnk';
elseif strcmp(opt.learning_type, 'edmd') 
    learning_type = 'edmd';
else
    learning_type = 'fro'; %default
end

% min and max of lambda range 
if isfield(opt, 'loglambda_range')
    loglambda_min = opt.loglambda_range(1);
    loglambda_max = opt.loglambda_range(2);
else   
    loglambda_min = 1e-8;
    loglambda_max = 1e10;
end

% min and max of rank range 
if isfield(opt, 'rnk')
    rnk_min = opt.rnk(1);
    rnk_max = opt.rnk(2);
else   
    rnk_min = 1;
    rnk_max = 20;
end

% number of Bayesian optimization iteration
if isfield(opt, 'nBO')
    nBO = opt.nBO;
else   
    nBO = 30;
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Bayesian optimization loop for estimating best hyperparameters, that are 
% either lambda or rank. This is done by CV.

% the case of Frobenius, operator, and nuclear norm regularization
if strcmp(learning_type, 'fro') || ...
          strcmp(learning_type, 'opr') || ...
            strcmp(learning_type, 'nuc') 
        
    rnk_var = optimizableVariable('loglambda', [loglambda_min,loglambda_max]);
    obj = @(theta)CV_LKO(theta,ZGY,learning_type);

    nBO_pp = 0;
    if 0
        results = bayesopt(obj, rnk_var, ...
            'IsObjectiveDeterministic',false,'Verbose', 0, ...
            'AcquisitionFunctionName','expected-improvement', ...
            'PlotFcn',{},'NumCoupledConstraints',0,...
            'MaxObjectiveEvaluations',nBO + nBO_pp);
    else
        results = bayesopt(obj,rnk_var, ...
            'IsObjectiveDeterministic',false,'Verbose',0,...
            'AcquisitionFunctionName','lower-confidence-bound',...
            'PlotFcn',{},'NumCoupledConstraints',0,...
            'MaxObjectiveEvaluations',nBO + nBO_pp);
    end

    if 0
        loglambda_est = results.XAtMinEstimatedObjective.loglambda;
    else
        loglambda_est = results.XAtMinObjective.loglambda;
    end
    %----------------------------------------------------------------------
% the case of rank constraint
elseif  strcmp(learning_type, 'rnk')
    rnk_max = min([rnk_max size(ZGY.KGtr,2)]);
    
    rnk_var = optimizableVariable('rnk', [rnk_min,rnk_max],'Type','integer');
    obj = @(theta)CV_LKO(theta,ZGY,learning_type);

    nBO_pp = 0;
    if 0
        results = bayesopt(obj, rnk_var, ...
            'IsObjectiveDeterministic',false,'Verbose', 0, ...
            'AcquisitionFunctionName','expected-improvement', ...
            'PlotFcn',{},'NumCoupledConstraints',0,...
            'MaxObjectiveEvaluations',nBO + nBO_pp);
    else
        results = bayesopt(obj,rnk_var, ...
            'IsObjectiveDeterministic',false,'Verbose',0,...
            'AcquisitionFunctionName','lower-confidence-bound',...
            'PlotFcn',{},'NumCoupledConstraints',0,...
            'MaxObjectiveEvaluations',nBO + nBO_pp);
    end

    if 0
        rnk_est = results.XAtMinEstimatedObjective.rnk;
    else
        rnk_est = results.XAtMinObjective.rnk;
    end
    %----------------------------------------------------------------------
% the case of EDMD    
else
    % lambda_est = 0;
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Main optimization for the learning problem

% the case of Frobenius norm regularization or rank constraint
if strcmp(learning_type, 'fro') || strcmp(learning_type, 'rnk')
    if strcmp(learning_type, 'fro') 
        lambda_est = exp(loglambda_est) ...
                         /(size(ZGY.Xt,1) * size(ZGY.KGt,1))...
                         *(size(ZGY.X,1)  * size(ZGY.KGt,1));
        A = opt_sol_fro(ZGY.Z,ZGY.KG,ZGY.Y,lambda_est);
    else
        [A,~,~] = opt_sol_rnk(ZGY.Z,ZGY.KG,ZGY.Y,rnk_est);
    end
    %----------------------------------------------------------------------
% the case of operator and nuclear norm regularization    
elseif  strcmp(learning_type, 'opr') || strcmp(learning_type, 'nuc') 
    lambda_est = exp(loglambda_est) ...
                /(size(ZGY.Xt,1) * size(ZGY.KGt,1))...
                *(size(ZGY.X,1)  * size(ZGY.KGr,1));
    B0 = (eye(size(ZGY.Zr))/(ZGY.Zr)) * (ZGY.Y) * (eye(size(ZGY.KGr))/(ZGY.KGr));     
  
    if strcmp(learning_type, 'nuc') 
        [A,~,~,~,~,~,~] = opt_sol_nuc(ZGY.Zr,ZGY.KGr,ZGY.Y,lambda_est,B0);
    else
        [A,~,~,~,~,~,~] = opt_sol_opr(ZGY.Zr,ZGY.KGr,ZGY.Y,lambda_est,B0);
    end
    %----------------------------------------------------------------------
% the case of EDMD    
else 
    A = opt_sol_fro(ZGY.Z,ZGY.KG,ZGY.Y,0);
    %----------------------------------------------------------------------
end
C = (ZGY.PGinv_GV) * A;
end
