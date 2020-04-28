function [gm,au] = getGravitationalParameters(s)
G = 6.67408e-11;
%m = [sun-1    mercury-2  venus-3  earth-4   mars-5] 
m = [1.9891e30 330.2e21 4.8685e24 6.0477e24 641.85e21];
gm = G*m(s);
au = 149597870.7e3;
end
