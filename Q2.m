%2.1a
clc
clear all;
m_sun = 1.9891e30; %kg
m_earth = 6.0477e24; %kg
mu = m_earth/(m_earth+m_sun);
au = 1.4959e11; %m
com_sunref = m_earth/(m_earth + m_sun);   %distance between sun and com
com_earthref = m_sun/(m_earth + m_sun);   %distance between earth and com
x = [1.010075186335748, -1.000001266837936];
y = 0;
r1 = sqrt((mu + x).^2);
r2 = sqrt((1 - mu - x).^2);
for i = 1:2
    Uxx(i) = 1 - ((1 - mu)/(r1(i)^3)) - (mu/(r2(i)^3)) + (3*(1-mu)*((x(i)+mu)^2)/(r1(i)^5)) + (3*mu*((x(i)+mu-1)^2)/(r2(i)^5));
    Uxy(i) = (3*(1-mu)*y*((x(i)+mu)^2)/(r1(i)^5)) + (3*mu*y*((x(i)+mu-1)^2)/(r2(i)^5));
    Uyy(i) = 1 - ((1 - mu)/(r1(i)^3)) - (mu/(r2(i)^3)) + (3*(1 - mu)*(y^2)/(r1(i)^5)) + (3*mu*(y^2)/(r2(i)^5));
    A{i} = [0 0 1 0; 0 0 0 1; Uxx(i) Uxy(i) 0 2; Uxy(i) Uyy(i) -2 0];
    [V{i} E{i}] = eig(A{i});
end

%% 
%Q2.1(b)
%if determinant of A is zero then the eigen values are correct
for i = 1:length(E)
    I = eye(length(A{i}));
    j=1;
    while j <= 4
        determinant{i}(j) = det(A{i} - E{i}(j,j)*I); %determinant = |A - lambda*I|
        j=j+1;
    end
end

%%
%2.2
%manifold for L2
pert = 1e-5;
for k = 1:length(V)
if real(E{1}(k,k)) > 0
V_L2 = real(V{1}(:,k)); %only the ones with eigen values exceeding zero
end
if real(E{1}(k,k)) < 0
V_2_L2 = real(V{1}(:,k)); %only the ones with eigen values lesser than zero
end
if real(E{1}(k,k)) > 0
V_3_L2 = real(V{1}(:,k)); %unstable manifold
end
if real(E{1}(k,k)) < 0
V_4_L2 = real(V{1}(:,k)); %stable manifold
end
end
for k = 1:length(V)
if real(E{2}(k,k)) > 0
V_L3 = real(V{2}(:,k)); %only the ones with eigen values exceeding zero
end
if real(E{2}(k,k)) < 0
V_2_L3 = real(V{2}(:,k)); %only the ones with eigen values lesser than zero
end
end
x_i = 1.010075186335748;
x_i_L3 = -1.000001266837936;
y_i = 0;
x_dot = 0;
y_dot = 0;
% 4 combinations of initial values for 4 different manifolds obtained by changing sign of perturbation and stability of manifold 
X_initial = [x_i y_i x_dot y_dot] + [-(pert*(V_L2(1)/norm(V_L2(1)))),-(pert*(V_L2(2)/norm(V_L2(2)))),-(pert*(V_L2(3)/norm(V_L2(3)))),-(pert*(V_L2(4)/norm(V_L2(4))))];
X_initial_L3 = [x_i_L3 y_i x_dot y_dot] + [-(pert*(V_L3(1)/norm(V_L3(1)))),-(pert*(V_L3(2)/norm(V_L3(2)))),-(pert*(V_L3(3)/norm(V_L3(3)))),-(pert*(V_L3(4)/norm(V_L3(4))))];
%X_initial = [x_i y_i x_dot y_dot] + [-(pert*V_L2(1)),-(pert*V_L2(2)),-(pert*V_L2(3)),-(pert*V_L2(4))];
t_i = 0;
t_f = 10*pi;
options = odeset('AbsTol',1e-12);
[~,X_output] = ode45(@loc,[t_i t_f],X_initial,options);

x_man = X_output(:,1);
y_man = X_output(:,2);

x_i_2 = 1.010075186335748;
x_i_2_L3 = -1.000001266837936;
y_i_2 = 0;
x_dot_2 = 0;
y_dot_2 = 0;

X_2_initial = [x_i_2 y_i_2 x_dot_2 y_dot_2] - [-(pert*(V_2_L2(1)/norm(V_2_L2(1)))),-(pert*(V_2_L2(2)/norm(V_2_L2(2)))),-(pert*(V_2_L2(3)/norm(V_2_L2(3)))),-(pert*(V_2_L2(4)/norm(V_2_L2(4))))];
X_2_initial_L3 = [x_i_L3 y_i x_dot y_dot] - [-(pert*(V_L3(1)/norm(V_L3(1)))),-(pert*(V_L3(2)/norm(V_L3(2)))),-(pert*(V_L3(3)/norm(V_L3(3)))),-(pert*(V_L3(4)/norm(V_L3(4))))];
%X_initial = [x_i y_i x_dot y_dot] + [-(pert*V_L2(1)),-(pert*V_L2(2)),-(pert*V_L2(3)),-(pert*V_L2(4))];
t_i = 0;
t_f = 10*pi;
options = odeset('AbsTol',1e-12);
[t,X_2_output] = ode45(@loc,[t_i -t_f],X_2_initial,options);

x_man_2 = X_2_output(:,1);
y_man_2 = X_2_output(:,2);

x_i_3 = 1.010075186335748;
y_i_3 = 0;
x_dot_3 = 0;
y_dot_3 = 0;

X_3_initial = [x_i_3 y_i_3 x_dot_3 y_dot_3] - [-(pert*(V_3_L2(1)/norm(V_3_L2(1)))),-(pert*(V_3_L2(2)/norm(V_3_L2(2)))),-(pert*(V_3_L2(3)/norm(V_3_L2(3)))),-(pert*(V_3_L2(4)/norm(V_3_L2(4))))];
%X_initial = [x_i y_i x_dot y_dot] + [-(pert*V_L2(1)),-(pert*V_L2(2)),-(pert*V_L2(3)),-(pert*V_L2(4))];
t_i = 0;
t_f = 10*pi;
options = odeset('AbsTol',1e-12);
[t,X_3_output] = ode45(@loc,[t_i t_f],X_3_initial,options);

x_man_3 = X_3_output(:,1);
y_man_3 = X_3_output(:,2);

x_i_4 = 1.010075186335748;
y_i_4 = 0;
x_dot_4 = 0;
y_dot_4 = 0;

X_4_initial = [x_i_4 y_i_4 x_dot_4 y_dot_4 ] + [-(pert*(V_4_L2(1)/norm(V_4_L2(1)))),-(pert*(V_4_L2(2)/norm(V_4_L2(2)))),-(pert*(V_4_L2(3)/norm(V_4_L2(3)))),-(pert*(V_4_L2(4)/norm(V_4_L2(4))))];
%X_initial = [x_i y_i x_dot y_dot] + [-(pert*V_L2(1)),-(pert*V_L2(2)),-(pert*V_L2(3)),-(pert*V_L2(4))];
t_i = 0;
t_f = 10*pi;
options = odeset('AbsTol',1e-12);
[t,X_4_output] = ode45(@loc,[t_i -t_f],X_4_initial,options);

x_man_4 = X_4_output(:,1);
y_man_4 = X_4_output(:,2);

figure;   
plot(x_man,y_man,'r');    %unstable
hold on
plot(x_man_2,y_man_2,'b');    %stable
hold on
plot(x_man_3,y_man_3,'r');   %unstable
hold on
plot(x_man_4,y_man_4,'b');   %stable
hold on
plot(0,0,'-o',1,0,'-o');
axis tight
legend('Unstable Manifold', 'Stable Manifold')
title('Manifolds for L2 point')
ylabel('Non-dimensionalised distance along Y axis')
xlabel('Non-dimensionalised distance along X axis')
grid on

figure
plot(x_man_4,y_man_4,'b');
hold on
plot(x_man_3,y_man_3,'r');
grid on

legend('Stable Manifold','Unstable Manifold');




