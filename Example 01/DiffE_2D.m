%--------------------------------------------------------------------------
% Learnin Koopman Operator
% Ref: Representer Theorem for Learning Koopman Operators 
% Link: https://arxiv.org/abs/2208.01681
%
% Dynamical system used in Example 1 of reference above.
%
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
function y = DiffE_2D(x)
mu1 = 0.95;
mu2 = 0.75;
mu  = mu1^2 - mu2;

y = [ mu1 * x(1);
      mu2 * x(2) + mu * x(1)^2];
end