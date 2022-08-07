%--------------------------------------------------------------------------
% Learnin Koopman Operator
% Ref: Representer Theorem for Learning Koopman Operators 
% Link: https://arxiv.org/abs/2208.01681
%
% The code simulate the dynamics for the discretized version of the following
% PDE describing u(xi,t) as 
%       du/dt = du/dxi + 0.1 * d^2u/dxi^2,  
% which is used in Example 4 in ref above.
%
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
function y = DiffE_PDE(u,Amx)
% A linear PDE with advection and diffusion term:
%     du/dt = du/dxi + 0.1 d^2u/dxi^2
% 
% discrete-time version
%     u_{t+1} = u_{t} + Dt (du/dxi + 0.1 d^2u/dx^2)
%
% discrete-time and discrete-space version
%     u_{t+1,k} = u_{t,k} + Dt [ 1/Dxi (u_{t,k+1}-u_{t,k}) ...
%                     +  0.1/Dxi^2 (u_{t,k+1}-2*u_{t,k}+u_{t,k-1}) ]
%               = u_{t,k} + Dt/Dxi (u_{t,k+1}-u_{t,k}) ...
%                     +  0.1 Dt/Dxi^2 (u_{t,k+1}-2*u_{t,k}+u_{t,k-1})
%               
% -------------------------------------------------------------------------

% boundary conditions
u(:,1) = 0;
u(:,end) = 0;

% using first and second discrete derivatives in xi-domain, we can write
% the discretized PDE as x_+ = F(x). 
v = u';
y = Amx * v;
y = y';
end