function f = CR3BP_2D(t,X)
global mu
x=X(1) ; y=X(2);
xdot = X(3); ydot = X(4);
r1 = sqrt((x+mu)^2+y^2); 
r2 = sqrt((x-1+mu)^2+y^2);
[Ux,Uy] = getdU(X(1:2));
f = [xdot;...
     ydot;...
     2*ydot+Ux;...
     -2*xdot+Uy];
end