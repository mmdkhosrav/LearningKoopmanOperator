%--------------------------------------------------------------------------
% The code symmetrize the matrix G, makes sure that it is positive definite 
% and also, returns its square root.
%
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
function [G,Gr] = mxPD(G,eps_G)
% to save a backup for G
% G0 = G;

% symmetrizing G
G  = 0.5*(G + G');

% replacing G with G + eps_G * I to slightly shift its eigenvalues for 
% ensuring its positive definiteness.
G = G + eps_G * eye(size(G));

% Gr: matrix square root of G
Gr = sqrtm(G);

% to avoid imaginary parts which can happen due to numerical imprecisions
Gr = real(Gr);

% symmetrizing Gr
Gr = 0.5*(Gr + Gr');

% ensuring Gr is the root of G
G  = Gr * Gr;

% Diff_GG0 = norm(G-G0);
end



