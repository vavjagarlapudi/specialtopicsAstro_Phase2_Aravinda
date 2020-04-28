clc; clear all;


global idx1
global idx2
global mu


idx1 = 3;

idx2 = 1; 5



tol = 1e-12; %tolerance
N = 100; %max no of iteration
x_0 = [1 -1]; %initial guess x
y_0 = [0 0];
mu = muCalculator(idx1);

[L, error] = getLagrange(mu, x_0, tol, N);
%L has lagrange points
[~, r1, r2] = getCR3BP(L, 0, mu);
U = 0.5*L.^2 + (1-mu)./r1 + mu./r2;

shift = [0.011 -0.011]; %shifts in lagrange point
x = L-shift; y = 0;
[f, r1, r2] = getCR3BP(x, y, mu);
n_cap(:,1) = -f(:,1)/norm(f(:,1));
n_cap(:,2) = -f(:,2)/norm(f(:,2));
for i = 1:length(x)
    r1_cap(i) = r1(:,i)/norm(r1(:,i));
    beta(:,i) =  ((r1(i).^2)./(1-mu)).*(norm(f(:,i))./(r1_cap(i)'.*n_cap(:,i)).^2);
end
n_cap(:,1) = -f(:,1)/norm(f(:,1));
n_cap(:,2) = -f(:,2)/norm(f(:,2));

acc = beta.*(1-mu)./r1.^2 * (dot(r1/norm(r1),[1 -1]))^2 .*n_cap; 
x = L; y = 0; z = 0;
for i = 1:length(x)
    [eigenVectors{i},eigenValues{i}, A{i}] = getEigen(x(i), y, z, mu);
end
%eigenValues contain the eigen values for  L2, L3


for i = 1:length(eigenValues)
    I = eye(length(A{i}));
    j=1;
    while j <= 4
        determ{i}(j) = det(A{i} - eigenValues{i}(j,j)*I); %|A - lambda*I|
        j=j+1;
    end
end
%determ contains A - lambdaI values, close to zero

x  = [L;0 0];
tmin = 0; tmax = 10*pi;
ms_tol = 1e-5;
ode_options = odeset('Reltol',1e-12,'AbsTol',1e-12);
figure
for i = 1:length(x)
    %getEigenVector(eigenVectors, eigenValues, stability)
    %stability = 1 for stable, 0 for unstable
    unstableEigenVector{i} = getEigenVector(eigenVectors{i}, eigenValues{i}, 0);
    unstable_ICs{i} = [x(:,i);0;0]+ms_tol*unstableEigenVector{i}/norm(unstableEigenVector{i});
    stableEigenVector{i} = getEigenVector(eigenVectors{i}, eigenValues{i}, 1);
    stable_ICs{i} = [x(:,i);0;0]+ms_tol*stableEigenVector{i}/norm(stableEigenVector{i});
    [t_unstable_2d{j},manifolds_unstable_2d{i}] = ode45(@CR3BP_2D,[tmin tmax],unstable_ICs{i},ode_options);
    [t_stable_2d{j},manifolds_stable_2d{i}] = ode45(@CR3BP_2D,[tmin -tmax],stable_ICs{i},ode_options);
    figure(i)
    grid on
    hold on
    plot(manifolds_unstable_2d{i}(:,1),manifolds_unstable_2d{i}(:,2),...
        manifolds_stable_2d{i}(:,1),manifolds_stable_2d{i}(:,2))
    xlabel('x (ND)'); ylabel('y (ND)');
    scatter(x(1,i),0, 'filled');
    if i==1
        idx = 1;
    elseif i==2
        idx = 3;
    end
    title(sprintf('L%d', idx));
    legend('Unstable Manifold', 'Stable Manifold', sprintf('L%d', idx));
end


%% Q 2.3 (a)

delta = 0; %clock angle
alpha = -70:2.5:70; %cone angle
assigned_beta = 0.1; %from question
tmin = 0; tmax = 10*pi;
for i=1:length(alpha)
    for j=1:length(x)
        unstable_state = [unstable_ICs{j}(1:2); 0; unstable_ICs{j}(3:4); 0];
        stable_state = [stable_ICs{j}(1:2); 0; stable_ICs{j}(3:4); 0];
        [t_unstable_3d{j},manifolds_unstable_3d{j}] = ode45(@(t,X) CR3BP_3D(t,X,assigned_beta,delta,alpha(i)),...
            [tmin tmax],unstable_state,ode_options);
        [t_stable_3d{j},manifolds_stable_3d{j}] = ode45(@(t,X) CR3BP_3D(t,X,assigned_beta,delta,alpha(i)),...
            [tmin -tmax],stable_state,ode_options);
        manifolds{i,j} = {manifolds_unstable_3d{j},manifolds_stable_3d{j}};
        times{i,j} = {t_unstable_3d{j}, t_stable_3d{j}};
        %key - manifolds contains two cells [unstable] [stable]
        % --->
        %each contains 2 more cells, [L1] [L3] 
    end
    % Euclidean Norm L1 ---> L3 %
    p2pdist = pdist2(real(manifolds_unstable_3d{1}),real(manifolds_stable_3d{2}), 'euclidean');
    min_p2pdist(i) = min(min(p2pdist));
    atalpha(i) = alpha(i);
    [L1i(i), L3i(i)] = find(p2pdist == min_p2pdist(i));
end
min_D=min(min_p2pdist); %min for least alpha
alpha_idx = find(min_p2pdist == min_D);
opt_alpha = alpha(alpha_idx);
%%
%%% plot for alpha = opt_alpha
figure
plot3(manifolds{alpha_idx(1),1}{1,1}(1:L1i(alpha_idx(1)),1),...
    manifolds{alpha_idx(1),1}{1,1}(1:L1i(alpha_idx(1)),2),...
    manifolds{alpha_idx(1),1}{1,1}(1:L1i(alpha_idx(1)),3));
hold on
plot3(manifolds{alpha_idx(1),2}{1,2}(1:L3i(alpha_idx(1)),1),...
    manifolds{alpha_idx(1),2}{1,2}(1:L3i(alpha_idx(1)),2),...
    manifolds{alpha_idx(1),2}{1,2}(1:L3i(alpha_idx(1)),3));
hold on
scatter3(manifolds{alpha_idx(1),1}{1,1}(L1i(alpha_idx(1)),1),...
    manifolds{alpha_idx(1),1}{1,1}(L1i(alpha_idx(1)),2),...
    manifolds{alpha_idx(1),1}{1,1}(L1i(alpha_idx(1)),3), 'r', 'filled');
hold on
scatter3(manifolds{alpha_idx(1),2}{1,2}(L3i(alpha_idx(1)),1),...
    manifolds{alpha_idx(1),2}{1,2}(L3i(alpha_idx(1)),2),...
    manifolds{alpha_idx(1),2}{1,2}(L3i(alpha_idx(1)),3), 'g', 'filled');
hold on
scatter3(x(1,1),0,0, '*')
hold on
scatter3(x(1,2),0,0, '*')
grid on
xlabel('x (ND)'); ylabel('y (ND)'); zlabel('z (ND)');
legend('Unstable Manifold', 'Stable Manifold', 'L2-exit', 'L3-entry', 'L2', 'L3');
title('Minimum Euclidean distance for L2 - L3 transfer')
%% 
au = 149597870.7;
%%% Calculation of position and velocity errors
pos_unstpt = norm(manifolds{alpha_idx(2),1}{1,1}(L1i(alpha_idx(2)),1:3));
pos_stpt = norm(manifolds{alpha_idx(2),2}{1,2}(L3i(alpha_idx(2)),1:3));
poserr = abs(pos_unstpt - pos_stpt)*au;
vel_unstpt = norm(manifolds{alpha_idx(2),1}{1,1}(L1i(alpha_idx(2)),4:6));
vel_stpt = norm(manifolds{alpha_idx(2),2}{1,2}(L3i(alpha_idx(2)),4:6));
velerr = abs(vel_unstpt - vel_stpt)*au/(2*pi*365.25*86400);
%%
%%%transfer times 
time_unst_L1 =  times{alpha_idx(1),1}{1,1}(L1i(alpha_idx(1)),1)/(2*pi);
time_st_L3 =  times{alpha_idx(1),2}{1,2}(L3i(alpha_idx(1)),1)/(2*pi);


%% Data input generation for ANN
%input_ANN syntax
% [alpha] [t_u] [t_s] [J = || Xu_tu - Xs_ts ||]
for i=1:length(alpha)
    input_ANN(i,1) = alpha(i);
    input_ANN(i,2) = times{i,1}{1,1}(L1i(i),1);
    input_ANN(i,3) = times{i,2}{1,2}(L3i(i),1);
    input_ANN(i,4) = min_p2pdist(i);
end
writematrix(input_ANN,'ANN_WP4.csv')
