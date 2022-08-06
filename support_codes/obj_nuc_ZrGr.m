%--------------------------------------------------------------------------
% The code returns the function
%       J(B) = ||Z^0.5 * B * G^0.5 - Y||_F^2 + lambda * ||B||_*
% and its sub-gradient.
%
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
function [y,g] = obj_nuc_ZrGr(x,Y,Zr,Gr,lambda)
%-------------------------------------------------------------------------%
% Y:    The matrix with entries g_l(x(k-1)), for l=1,...,nG and k=1,...,nS  
% Zr:   Z^0.5 
% Gr:   G^0.5 
%--------------------------------------------------------------------------
eps_B = 1e-12;

nZ = size(Zr,1);
nG = size(Gr,1);
nZG = min([nZ,nG]);

B = reshape(x,nZ,nG);

E = Zr * B * Gr - Y;

dB = zeros(nZ,nG);
dB(1:nZG,1:nZG) = eps_B * eye(nZG); 
[~,S,~] = svd(B + dB);
y = sum(E(:).^2) + lambda * sum(diag(S));

%--------------------------------------------------------------------------
if nargout == 2
    eps_B = 1e-10;
    dB = zeros(nZ,nG);
    dB(1:nZG,1:nZG) = eps_B * eye(nZG); 
    [U,~,V] = svd(B + dB);
    Gv = 2 * Zr * E * Gr + lambda * U(:,1:nZG) * V(:,1:nZG)';
    g =   Gv(:);
end

end