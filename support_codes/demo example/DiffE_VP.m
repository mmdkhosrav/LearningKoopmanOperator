%--------------------------------------------------------------------------
% Learnin Koopman Operator
% Ref: Representer Theorem for Learning Koopman Operators 
% Link: https://arxiv.org/abs/2208.01681
%
% The Van der Pol dynamics used in Example 2 of the reference above.
%
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
function y = DiffE_VP(x,dt)
eps = 0.5;

dx1 = x(:,1) + dt * x(:,2);
dx2 = x(:,2) + dt *( eps * (1 - x(:,1).^2).*x(:,2) - x(:,1) );

y = [dx1 dx2];
end