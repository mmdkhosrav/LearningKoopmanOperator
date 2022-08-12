% -------------------------------------------------------------------------
% Nicholson and Bailey model for host and parasitoids is defined as
%       H_{n+1} = R0 H_n f(H_n, P_n),       
%       P_{n+1} = c H_n (1 - f (H_n, P_n)),
% where f is as follows
%       f(H,P) = (1+b P/k)^-k.
% 
% We have used the change-of-variable defined as follows
%       x1_n = alpha * H_n,
%       x2_n = beta  * P_n,
% where alpha = beta = 0.001.
%
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
function y = DiffE_HP(x)

R0 = 1.1; 
b  = 0.001; 
c  = 3;
k = 0.5;

alpha = 1e-3; 
beta =   1e-3;

f = (1 + b*x(:,2)/(beta * k)).^(-k);

y1 = alpha * R0 * x(:,1) .* f/alpha; 
y2 = beta  * c  * x(:,1) .* (1 - f)/beta; 

y = [y1 y2];
end