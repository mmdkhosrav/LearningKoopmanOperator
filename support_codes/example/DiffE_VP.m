function y = DiffE_VP(x,dt)
% mu = 0.8;
% a = 0.7; 
% b = 0.8; 
% c = 0.08; 
% 
% dx1 = x(:,1) + dt * ( x(:,1) +  x(:,1).^3/3 - x(:,2) + mu);
% dx2 = x(:,2) + dt * c * (x(:,2) + x(:,1)- b*x(:,2) + a );


eps = 0.5;

dx1 = x(:,1) + dt * x(:,2);
dx2 = x(:,2) + dt *( eps * (1 - x(:,1).^2).*x(:,2) - x(:,1) );

y = [dx1 dx2];
end