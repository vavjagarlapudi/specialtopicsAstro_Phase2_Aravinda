function real_V = getEigenVector(eigenVectors, eigenValues, stability)
%function that extracts alpha < 0 (stable/unstable lagrange points)
if stability == 1
    check = real(eigenValues) < 0 ;
elseif stability ==0
    check = real(eigenValues) > 0 ;
end
for i = 1:length(check)
    if check(i,i) == 1
        idx = i;
    end
end
real_V = eigenVectors(:,idx);
end