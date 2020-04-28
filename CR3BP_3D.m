function f = CR3BP_3D(t,X,beta,delta,alpha)
global idx1
global mu
x=X(1); 
y=X(2); 
z=X(3);
x_dot = X(4); 
y_dot = X(5); 
z_dot = X(6);
r1 = sqrt((x+mu)^2+y^2+z^2); 
r2 = sqrt((x-1+mu)^2+y^2+z^2);
r1_vector = [mu + x; y; z]; 
Ux = x - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;
Uy = y - (1-mu)*y/r1^3 - mu*y/r2^3;
Uz = -(1-mu)*z/r1^3-mu*z/r2^3;
f1 = [x_dot;...
        y_dot;...
        z_dot;...
        2*y_dot+Ux;...
        -2*x_dot+Uy;...
        Uz];
n = getnvector(alpha, delta);
[r1_cap, theta_cap, eta_cap] = getSolarSailRefFrame(r1_vector);
n_cap = [r1_cap theta_cap eta_cap]*n;
a_solar_sail = beta*(1-mu)/r1^2 * (r1_cap'*n_cap)^2 * n_cap;
f = f1 + [zeros(3,1); a_solar_sail];
end