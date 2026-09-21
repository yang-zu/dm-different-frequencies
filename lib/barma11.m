function [u,v] = barma11(phi1,theta1,phi2,theta2,rho,n)
% simulate error process following bivariate arma(1,1) model
% mean is 0 always

eps = randn(n+1,2);
eps(:,2)= rho*eps(:,1)+sqrt(1-rho^2)*eps(:,2);

u = zeros(n,1);
v = zeros(n,1);
u(1) = eps(1,1);
v(1) = eps(1,2);

for i=2:n
    u(i) = phi1*u(i-1)+eps(i,1)+theta1*eps(i-1,1);
    v(i) = phi2*v(i-1)+eps(i,2)+theta2*eps(i-1,2);
end



end