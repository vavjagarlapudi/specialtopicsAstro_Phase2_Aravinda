function mu = muCalculator(b)
m_sun = 1.9891e30;
%m = [mercury-1  venus-2  earth-3   mars-4] 
m = [330.2e21 4.8685e24 6.0477e24 641.85e21];
mu = m(b)/(m(b) + m_sun);
end