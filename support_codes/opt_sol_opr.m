%--------------------------------------------------------------------------
% The code solves optimization problem
%      min_A || Z * A * G - Y ||_F^2 + lambda * || Z^0.5 * A * G^0.5 ||^2
% The solution can be obtained using L-BFGS.
%
% Requirements: obj_opr_ZrGr.m 
% Potential requirements: fminlbfgs.m
% 
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
function [AB_star, B_star , fvalB, time_elapsed, exitflag, output, grad] =...
    opt_sol_opr(Zr,Gr,Y,lambda,B0)
%-------------------------------------------------------------------------%
% Y:    The matrix with entries g_l(x(k-1)), for l=1,...,nG and k=1,...,nS  
% Zr:   Z^0.5 
% Gr:   G^0.5 
%--------------------------------------------------------------------------
nZ = size(Zr,1);
nG = size(Gr,1);

if (nargin==4)
   B0 = (eye(nZ)/Zr) * Y * (eye(nG)/Gr);
end
x0 = B0(:);

% determining the solver choice, 0: fmincon with BFGS default, and 1: L-BFGS 
% borrowed from the following link:
% www.mathworks.com/matlabcentral/fileexchange/23245-fminlbfgs-fast-limited-memory-optimizer
solver_choice = 1;

if solver_choice
    tic
    opt_LBFGS = struct(  'GradObj' , 'on' ...
                     , 'Display' , 'off' ...'off''iter'
                     , 'LargeScale' , 'on' ...
                     , 'HessUpdate' , 'lbfgs'...
                     , 'StoreN' , 15 ...
                     , 'InitialHessType' , 'identity'...
                     , 'GoalsExactAchieve' , 1 ...
                     , 'GradConstr' , false ...
                     , 'TolX' , 1e-8 ...
                     , 'TolFun' , 1e-8 ...
                     , 'MaxIter' , 1e8 ...
                     , 'MaxFunEvals', 1e8 ...
                 );

    [x,fvalB,exitflag,output,grad] = ...
        fminlbfgs(@(x)obj_opr_ZrGr(x,Y,Zr,Gr,lambda),x0,opt_LBFGS);
else
    opt_fminunc = optimoptions('fminunc'...
                     , 'Display' , 'off' ...'off''iter'
                     , 'SpecifyObjectiveGradient',true);

    [x,fvalB,exitflag,output,grad] = ...
        fminunc(@(x)obj_opr_ZrGr(x,Y,Zr,Gr,lambda),x0,opt_fminunc);
end

B_star = reshape(x,nZ,nG);

AB_star = (eye(nZ)/Zr) * B_star * (eye(nG)/Gr);
time_elapsed = toc;
end