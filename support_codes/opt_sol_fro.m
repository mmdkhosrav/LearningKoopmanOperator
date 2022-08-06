%--------------------------------------------------------------------------
% The code solves optimization problem
%      min_A || Z * A * G - Y ||_F^2 + lambda * || Z^0.5 * A * G^0.5||_F^2,
% which has a closed form solution.
% 
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
function A_star = opt_sol_fro(Z,G,Y,lambda)
%--------------------------------------------------------------------------
eps_sol = 1e-12;

[U_Z, S_Z, V_Z] = svd(Z);
U = 0.5*(U_Z+V_Z);  
clear U_Z V_Z           % to save space on RAM
s_Z = diag(S_Z);    
clear S_Z               % to save space on RAM

[U_G, S_G, V_G] = svd(G);
V = 0.5*(U_G+V_G); 
clear U_G V_G           % to save space on RAM
s_G = diag(S_G);   
clear S_G               % to save space on RAM

W = U' * Y * V;
A_star = U * (W./(lambda + eps_sol + s_Z*s_G')) * V';
end