%--------------------------------------------------------------------------
% The code returns the function
%       J(B) = ||Z^0.5 * B * G^0.5 - Y||_F^2 + lambda * ||B||^2
% and its sub-gradient.
%
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
function [y,g] = obj_opr_ZrGr(x,Y,Zr,Gr,lambda)
%--------------------------------------------------------------------------
% Y:    The matrix with entries g_l(x(k-1)), for l=1,...,nG and k=1,...,nS  
% Zr:   Z^0.5 
% Gr:   G^0.5 
%--------------------------------------------------------------------------
eps_B = 1e-12;

nZ = size(Zr,1);
nG = size(Gr,1);
nZG = min([nZ,nG]);

B = reshape(x,nZ,nG);

E = Zr  * B * Gr - Y;

dB = zeros(nZ,nG);
dB(1:nZG,1:nZG) = 0 * eps_B * eye(nZG); 
if nZG < 150
    [~,S,~] = svd(B + dB);
else
    [~,S,~] = svds(B + dB,1);   
end
y = sum(E(:).^2) + lambda * S(1)^2;

%--------------------------------------------------------------------------
if nargout == 2
    eps_B = 1e-12;
    dB = zeros(nZ,nG);
    dB(1:nZG,1:nZG) = eps_B * eye(nZG); 
    if nZG < 2500
        [U,S,V] = svd(B + dB);
    else
        [U,S,V] = svds(B + dB,1);
    end
    
    dGv = 2 * lambda * S(1) * U(:,1) * V(:,1)';
    s = diag(S);
    eps_eig = 1e-8;
    idx_s = find((s>(s(1) - eps_eig)));
    for k = idx_s(2:end)
        dGv = dGv + 2 * lambda * s(1) * U(:,k) * V(:,k)';
    end    
    Gv = 2 * Zr * E * Gr + dGv;

    g =   Gv(:);
end

end