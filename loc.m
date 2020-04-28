function F = loc(t,X)
x = X(1);
y = X(2);
u = X(3);
v = X(4);
m_sun = 1.9891e30; %kg
m_earth = 6.0477e24; %kg
mu = m_earth/(m_earth+m_sun);
r1 = sqrt((mu + x)^2+y^2);
r2 = sqrt((1 - mu - x)^2+y^2);
p1 = (1 - mu) .* (mu + x);
p2 = mu .* (1-mu-x);
F= [u; v;...
    2*v+x - p1./r1.^3 + p2./r2.^3;...
    -2*u+y - (1-mu).*y./r1.^3 - mu.*y./r2.^3];
end
