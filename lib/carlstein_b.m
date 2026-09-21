function m=carlstein_b(phi,n)
% this function compute the Carlstein block size given the AR coefficient
% to be phi

a = 1-phi;
b = 1+phi;
c = a*b;
m = (((2*phi/c)^2)^(1/3))*n^(1/3);
m = ceil(m);



end