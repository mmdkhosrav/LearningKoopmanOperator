%--------------------------------------------------------------------------
% The code returns the function
%       J(B) = ||Z^0.5 * B * G^0.5 - Y||_F^2 + lambda * ||B||_F^2
% and its gradient.
%
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
function [y,g] = obj_fro_ZrGr(x,Y,Zr,Gr,lambda)
%--------------------------------------------------------------------------
% Y:    The matrix with entries g_l(x(k-1)), for l=1,...,nG and k=1,...,nS  
% Zr:   Z^0.5 
% Gr:   G^0.5 
%--------------------------------------------------------------------------
nZ = size(Zr,1);
nG = size(Gr,1);

A = reshape(x,nZ,nG);

E = Zr  * A * Gr - Y;

y = sum(E(:).^2) + lambda * sum(A(:).^2);

%--------------------------------------------------------------------------
if nargout == 2
    
    Gv = 2 * Zr * E * Gr + 2 * lambda * A;
    g =  Gv(:);
end

end