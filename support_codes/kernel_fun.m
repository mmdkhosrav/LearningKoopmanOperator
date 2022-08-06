%--------------------------------------------------------------------------
% The code returns inner products of sections of kernel at the rows of x1 
% and x2. More precisely, if x1 = [x_11' ; x_12' ; ... ; x_1n1']' and  x2 = 
% [x_21' ; x_22' ; ...;x_2n2']', where x_1i and x_2j belongs to R^d, for i 
% 1, 2,..., n1 and j = 1, 2,..., n2, then Kx1x2 is a n1 by n2 matrix such 
% that  Kx1x2_ij = k(x_1i,x_2j).
% 
% Mohammad Khosravi
% Email: mohammad.khosravi@tudelft.nl
% Delft Center for Systems and Control (DCSC)
% Delft University of Technology (TU Delft) 
% August 2022
%--------------------------------------------------------------------------
function Kx1x2 = kernel_fun(x1,x2,theta,kernel_type)
% x1:           input location variable 1
% x2:           input location variable 2
% theta:        hyperparameters
% kernel_type:  kernel types: 
%                   SE: Squared Exponential or Gaussian kernel
%                   RQ: Rational Quadratic kernel
%                   P:  Periodic kernel
%                   LP: Locally Periodic kernel
%                   R:  Radial kernel
%                   M3: Matern kernel 3/2
%                   M5: Matern kernel 5/2
%                   L:  Linear kernel
%                   PN: Polynomial kernel
%
% Kx1x2:        The matrix of inner products of sections of kernel at the 
%               rows of x1 and x2. 
%               dim(K_x) = dim(x1) x dim(x2)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
nx1 = size(x1,1); % # of input 1
nx2 = size(x2,1); % # of input 2
nd  = size(x1,2); % dimension of input


if strcmp(kernel_type,'SE')
    P = zeros(nx1,nx2);
    for i=1:nx1
       P(i,:) = sum(((x2 - x1(i,:)).^2),2); 
    end
    l   = theta(1);
    Kx1x2 = exp(-P/(2*l^2));    
    
elseif strcmp(kernel_type,'RQ')
    P = zeros(nx1,nx2);
    for i=1:nx1
       P(i,:) = sum(((x2 - x1(i,:)).^2),2); 
    end
    l     = theta(1);
    alpha = theta(2);
    Kx1x2   = (1 + P/(2*alpha*l^2)).^(-alpha);    

elseif strcmp(kernel_type,'P')
    l   = theta(1);
    p   = theta(2);    
    D = ones(nx1,nx2);
    for j = 1:nd
        P = zeros(nx1,nx2);
        for i=1:nx1
           P(i,:) = abs( x2(:,j) - x1(i,j)); 
        end
        P = exp(-2*(sin(pi*P/p)).^2 ./(l^2));
        D = D .* P;
    end    
    Kx1x2 = D;    
    
elseif strcmp(kernel_type,'LP')
    l   = theta(1);
    p   = theta(2);    
    P = zeros(nx1,nx2);
    for i=1:nx1
       P(i,:) = sum(((x2 - x1(i,:)).^2),2); 
    end
    D = exp(-P/(2*l^2));    
    for j = 1:nd
        P = zeros(nx1,nx2);
        for i=1:nx1
           P(i,:) = abs( x2(:,j) - x1(i,j)); 
        end
        P = exp(-2*(sin(pi*P/p)).^2 ./(l^2));
        D = D .* P;
    end    
    Kx1x2 = D;    

elseif strcmp(kernel_type,'R')
    P = zeros(nx1,nx2);
    for i=1:nx1
       P(i,:) = sum(((x2 - x1(i,:)).^2),2); 
    end
    R   = P.^0.5;
    l   = theta(1);
    Kx1x2 = exp(-R/l);

elseif strcmp(kernel_type,'M3')
    P = zeros(nx1,nx2);
    for i=1:nx1
       P(i,:) = sum(((x2 - x1(i,:)).^2),2); 
    end
    R   = P.^0.5;
    l   = theta(1);
    Kx1x2 = (1 + 3^0.5*R/l) .* exp(-3^0.5*R/l);

elseif strcmp(kernel_type,'M5')
    P = zeros(nx1,nx2);
    for i=1:nx1
       P(i,:) = sum(((x2 - x1(i,:)).^2),2); 
    end
    R   = P.^0.5;
    l   = theta(1);
    Kx1x2 = (1 + 5^0.5*R/l + 5*P/(3*l^2)) .* exp(-5^0.5*R/l);

elseif strcmp(kernel_type,'PN')
    sigma_b = theta(1);
    sigma_p = theta(2);
    c = theta(end-nd+1:end);
    c = reshape(c,1,[]);    
    
    P = zeros(nx1,nx2);
    for i=1:nx2
       P(:,i) = (x1 - c)* (x2(i,:)-c)'; 
    end
    Kx1x2  = (sigma_b^2 + P)^sigma_p;  

elseif strcmp(kernel_type,'L')
    sigma_b = theta(1);
    c = theta(end-nd+1:end);
    c = reshape(c,1,[]);    
    
    P = zeros(nx1,nx2);
    for i=1:nx2
       P(:,i) = (x1 - c)* (x2(i,:)-c)'; 
    end
    Kx1x2  = sigma_b^2 + P;  
    
elseif strcmp(kernel_type,'LSE')
    P = zeros(nx1,nx2);
    for i=1:nx1
       P(i,:) = sum(((x2 - x1(i,:)).^2),2); 
    end
    l   = theta(1);

    K_x_SE = exp(-P/(2*l^2));
    
    sigma_b = theta(2);
    c = theta(end-nd+1:end);
    c = reshape(c,1,[]);    
    P = zeros(nx1,nx2);
    for i=1:nx2
       P(:,i) = (x1 - c)* (x2(i,:)-c)'; 
    end
    K_x_L  = sigma_b^2 + P;
    
    Kx1x2 = K_x_SE .* K_x_L; 
    
else
    sigma_b = theta(1);
    P = zeros(nx1,nx2);
    for i=1:nx2
       P(:,i) = x1*x2(i,:)'; 
    end
    Kx1x2  = sigma_b^2 + P;  

end

end