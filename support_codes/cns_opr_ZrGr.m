%--------------------------------------------------------------------------
% The code returns the constraint
%       norm(B) <= (1-eps_stb)  
% and its sub-gradient.
%
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
function [c,ceq,gradc,gradceq]  = cns_opr_ZrGr(x,nZ,nG,eps_stb)

nZG = min([nZ,nG]);

A = reshape(x,nZ,nG);

if nZG < 150
    [~,S,~] = svd(A);
else
    [~,S,~] = svds(A,1);    
end
c = S(1)^2 - (1 - eps_stb)^2;

ceq = [];
% -------------------------------------------------------------------------
if nargout > 2
    
    eps_A = 1e-11;
    dA = zeros(nZ,nG);
    dA(1:nZG,1:nZG) = eps_A * eye(nZG); 
    if nZG < 150
        [U,S,V] = svd(A + dA);
    else
        [U,S,V] = svds(A + dA,1);
    end
    
    Gradc = 2 * S(1) * U(:,1) * V(:,1)';
    s = diag(S);
    eps_eig = 1e-8;
    idx_s = find((s>(s(1) - eps_eig)));
    for k = idx_s(2:end)
        Gradc = Gradc + 2 * s(1) * U(:,k) * V(:,k)';
    end    
    gradc =  Gradc(:);    
    
    gradceq = [];
end

end