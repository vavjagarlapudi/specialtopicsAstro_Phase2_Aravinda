function [f,r1,r2] = getCR3BP(x, y, mu)
r1 = sqrt((mu + x).^2);
r2 = sqrt((1 - mu - x).^2);
p1 = (1 - mu) .* (mu + x);
p2 = mu .* (1-mu-x);
f = [x - p1./r1.^3 + p2./r2.^3;...
    y - (1-mu).*y./r1.^3 - mu.*y./r2.^3];
end
