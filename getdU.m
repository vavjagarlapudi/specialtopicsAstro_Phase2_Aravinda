function [Ux,Uy] = getdU(X)
global idx1
mu = muCalculator(idx1);
x=X(1) ; y=X(2);
r1 = sqrt((x+mu)^2+y^2); r2 = sqrt((x-1+mu)^2+y^2);
Ux = x - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;
Uy = y - (1-mu)*y/r1^3 - mu*y/r2^3;
end