%computation of required sail surface area
%1.2d
clc
clear all
m = 10; %mass of spacecraft in kg
L = 3.828e26;   %Luminosity in W
alpha = 0;    %angle of sail wrt radiation vector
au = 1.4959e11; %m
%distance from sun in m
R = [0.999075186335748+3.040411046456244e-06 -0.989001266837936+3.040411046456244e-06]*au;
c = 2.99792458e8; %speed of light in m/s
mu_sun = 1.3271e20;   %m^3/s^2
beta = [3.569014802472770 0.032636868519920]; %change if original script is changed
for i = 1:2
    P(i) = 2 * L * (cosd(alpha)^2)/(4 * pi * (R(i)^2) * c);   %radiation pressure in Pa
    sigma_lim = 1.53;   %
    A(i) = 1000 * m * beta(i)/sigma_lim;
    side(i) = sqrt(A(i));
end


