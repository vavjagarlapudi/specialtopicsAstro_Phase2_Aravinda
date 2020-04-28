function [eigenVectors, eigenValues, A] = getEigen(x, y, z, mu)
[~, r1, r2] = getCR3BP(x, y, mu);
mu1 = 1 - mu;
mu2 = mu;
Abar = mu1./r1.^3 + mu2./r2.^3;
Bbar = 3* (mu1./r1.^5 + mu2./r2.^5);
Cbar = 3*(mu1.*(x + mu2)./r1.^5 + mu2.*(x - mu1)./r2.^5);
Dbar = 3*(mu1.*(x + mu2).^2./r1.^5 + mu2.*(x - mu1).^2./r2.^5);
Uxx = 1 - Abar + Dbar;
Uyy = 1 - Abar + Bbar*y^2;
Uxy = Cbar*y;
Uzz = -Abar + Bbar*z^2;
Uxz = Cbar*z;
Uyz = Bbar*y*z;
A = [0 0 1 0;0 0 0 1; Uxx Uxy 0 2; Uxy Uyy -2 0];
[eigenVectors, eigenValues] = eig(A);
end
