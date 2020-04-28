function [L, error] = getLagrange(mu, x_0, tol, N)
%% CR3BP eqn
syms x y
f = getCR3BP(x, y, mu); %co ordinates needs to be in sym

%% NR Method for collinear Lagrange Points
x0 = x_0;
g = diff(f(1));
for j = 1:length(x0)
    for i=1:N
        f0=vpa(subs(f(1),x,x0(j)));
        df0=vpa(subs(g,x,x0(j)));
        val=x0(j)-f0/df0;
        err(i)=val-x0(j);
        if abs(err(i)) < tol
            break
        end
        x0(j)=val;
    end
    L(j) = double(val);
    error{j} = err;
end
end

