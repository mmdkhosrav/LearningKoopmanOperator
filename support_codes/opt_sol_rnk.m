%--------------------------------------------------------------------------
% The code solves optimization problem
%      min_A    || Z * A * G - Y ||_F^2 
%       s.t.    rank( Z^0.5 * A * G^0.5 ) <= r,
% which has a closed form solution.
% 
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
function [A_star, S, time_elapsed] = opt_sol_rnk(Z,G,Y,r)
%--------------------------------------------------------------------------
tic

nZ = size(Z,1);
nG = size(G,1);
nZG = min([nZ nG]);

r = floor(r);
r = min([r nZG]);

[U,S,V] = svd(Y);
S_tr = S;

if max([nZ nG]) < 500
    s = diag(S);
    s_tr = s;
    s_tr(r+1:end) = 0;
    S_tr(1:nZG,1:nZG) = diag(s_tr);
    Y_tr = U * S_tr * V';

    eps_A = 1e-6;
    A_star = (Z + eps_A * eye(nZ)) \ ...
            (Y_tr * (eye(nG) / (G + eps_A * eye(nG))));
else
    s = diag(S);
    s_tr = s;
    s_tr(r+1:end) = 0;
    S_tr(1:nZG,1:nZG) = diag(s_tr);
    eps_A = 1e-6;
    A_star = (Z + eps_A * eye(nZ)) \ ...
            ((U * S_tr * V') * (eye(nG) / (G + eps_A * eye(nG))));    
end

time_elapsed = toc;
end