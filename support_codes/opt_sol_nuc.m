%--------------------------------------------------------------------------
% The code solves optimization problem
%      min_A || Z * A * G - Y ||_F^2 + lambda * || Z^0.5 * A * G^0.5 ||_*
% The solution can be obtained using L-BFGS.
%
% Requirements: obj_nuc_ZrGr.m 
% Potential requirements: fminlbfgs.m
% 
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
function [AB_star, B_star , fvalB, time_elapsed, exitflag, output, grad] =...
    opt_sol_nuc(Zr,Gr,Y,lambda,B0)
%-------------------------------------------------------------------------%
% Y:    The matrix with entries g_l(x(k-1)), for l=1,...,nG and k=1,...,nS  
% Zr:   Z^0.5 
% Gr:   G^0.5 
%--------------------------------------------------------------------------
nZ = size(Zr,1);
nG = size(Gr,1);

if (nargin==4) || isempty(B0)
   B0 = (eye(nZ)/Zr) * Y * (eye(nG)/Gr);
end
x0 = B0(:);

tic
opt_LBFG = struct(  'GradObj' , 'on' ...
                 , 'Display' , 'off' ...'off''iter'
                 , 'LargeScale' , 'on' ...
                 , 'HessUpdate' , 'lbfgs'...
                 , 'StoreN' , 1 ...
                 , 'InitialHessType' , 'identity'...
                 , 'GoalsExactAchieve' , 1 ...
                 , 'GradConstr' , false ...
                 , 'TolX' , 1e-12 ...
                 , 'TolFun' , 1e-12 ...
                 , 'MaxIter' , 1e12 ...
                 , 'MaxFunEvals', 1e12 ...
             );
[x,fvalB,exitflag,output,grad] = fminlbfgs(@(x)obj_nuc_ZrGr(x,Y,Zr,Gr,lambda),x0,opt_LBFG);       

% opt_fminunc = optimoptions('fminunc'...
%                  , 'Display' , 'off' ...'off''iter'
%                  , 'SpecifyObjectiveGradient',true);
%          
% [x,fvalB,exitflag,output,grad] = fminunc(@(x)obj_nuc_ZrGr(x,Y,Zr,Gr,lambda),x0,opt_fminunc);

B_star = reshape(x,nZ,nG);

AB_star = (eye(nZ)/Zr) * B_star * (eye(nG)/Gr);
time_elapsed = toc;
end