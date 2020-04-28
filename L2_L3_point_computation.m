%%
%1a
%computation of collinear lagrange points L2 and L3 for Earth Sun 
%y = 0
%dU/dx = 0
clc
clear all;
m_sun = 1.9891e30; %kg
m_earth = 6.0477e24; %kg
mu = m_earth/(m_earth+m_sun);
au = 1.4959e11; %m
com_sunref = m_earth/(m_earth + m_sun);   %distance between sun and com
com_earthref = m_sun/(m_earth + m_sun);   %distance between earth and com
c = 1;
for j = -1.5:0.1:1.5   %initial guess
    x(1) = j;
for i = 1:1000
    r1_norm(i) = sqrt((mu + x(i))^2);   %r1 is from the sun
    r2_norm(i) = sqrt((1 - mu - x(i))^2);   %r2 is from the earth
    f_x(i) = x(i) - (((1-mu)/(r1_norm(i)^3))*(x(i) + mu)) - ((mu/(r2_norm(i)^3))*(x(i) - (1 - mu)));
    diff_f(1) = 1;
    h(1) = 10;
    if i>1
        h(i) = x(i) - x(i-1);
        diff_f(i) = (f_x(i)-f_x(i-1))/h(i);   %df(x)/dx
        if abs(x(i) - x(i-1)) < 1e-8
            x_l(c) = x(i);
            x_initial(c) = j; 
            r1_norm_L(c) = com_sunref + x_l(c);
            r2_norm_L(c) = -com_earthref + x_l(c);
            c = c + 1;  
            break
        end
    end    
    if i < 1000
    x(i+1) = x(i) - f_x(i)/diff_f(i);
    end
  
end
sol = [ x_initial' x_l' r1_norm_L' r2_norm_L'];
end

%%
%computation of required lightness number for shifted L2 and L3
for k = 1:length(x_l)    %we know L3 is behind sun and barycenter along x
    if x_l(k) > 1    %we know L2 is after earth
        L2 = x_l(k);
        L2_shift = L2 - 0.011;
        r1_L2 = r1_norm_L(k);
        r2_L2 = r2_norm_L(k);
    end
    
    if x_l(k) < 0
        L3 = x_l(k);
        L3_shift = L3 + 0.011;
        r1_L3 = r1_norm_L(k);
        r2_L3 = r2_norm_L(k);
    end
end

%computation of lightness number for L2
r1_L2_art = mu + L2_shift;   %r1 for the artificial L2
r2_L2_art = 1 - mu - L2_shift;   %r2 for the artificial L2
r1_L2_art_cap = r1_L2_art/norm(r1_L2_art); %unit vector for r1_L2_art
%dU/dx only because the rest are 0
grad_U_L2 = sqrt(L2_shift^2) - ((1 - mu)*(mu+L2_shift)/((r1_L2_art)^3)) + (mu*(1-mu-L2_shift)/((norm(r2_L2_art)^3)));
%orientation vector for the sail- along the direction of the sun only
n_L2 = grad_U_L2/norm(grad_U_L2);   
%solar sail lightness number
beta_L2 = (norm(r1_L2_art)^2)*grad_U_L2*n_L2/((1 - mu)*(norm(r1_L2_art_cap)*n_L2)^2);

%computation of lightness number for L3
L3_12 = [L3_shift 0 0];
r1_L3_art = mu + L3_shift;
r2_L3_art = [1 - mu - L3_shift 0 0];
r1_L3_art_cap = (r1_L3_art/norm(r1_L3_art));
grad_U_L3 = sqrt(L3_shift^2) - ((1 - mu)*(mu+L3_shift)/(r1_L3_art)^3) + (mu*(1-mu-L3_shift)/(norm(r2_L3_art)^3));
n_L3 = grad_U_L3/norm(grad_U_L3);
beta_L3 = (norm(r1_L3_art)^2)*grad_U_L3*n_L3/((1 - mu)*(norm(r1_L3_art_cap)*n_L3)^2);

%%
%solar sail acceleration for verification
beta = [beta_L2 beta_L3];
n = [n_L2 n_L3];
r1_art = [r1_L2_art r1_L3_art];
grad_U = [grad_U_L2 grad_U_L3];
r1_cap = [r1_L2_art_cap,r1_L3_art_cap];
for i = 1:2
a_sail(i) = beta(i)* (1 - mu)* (dot(r1_cap(i), n(i))^2)* n(i)/(r1_art(i)^2);
ver = a_sail - grad_U;
end