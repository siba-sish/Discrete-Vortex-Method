% Code for flow physics analysis of a cambered airfoil using the Discrete Vortex Method
clc;
clear;

syms f(x)
y = f(x);
y = 4*0.05*x*(1-x);

%Variable declaration and initialisation
k= 0.05*(4-8*x);
N=10; % Number of divisions
c=1;
eps = 0.05;
x = linspace(0,c,N+1);
x_01 = zeros(N,1);
x_c1 = zeros(N,1);
z_01 = zeros(N,1);
z_c1 = zeros(N,1);
a = zeros(N);
dx = x(2) - x(1);

x_01(1) = 0.25*dx;
x_c1(1) = 0.75*dx;

% Coordinates of required points in panels
for  i =2:N;
   x_01(i) = dx + x_01(i-1);
   x_c1(i) = dx + x_c1(i-1);
   z(i) = 4*x(i)*eps*(1-x(i));
   z_01(i) = 4*x_01(i)*eps*(1-x_01(i));
   z_c1(i) = 4*x_c1(i)*eps*(1-x_c1(i));
end
z(N+1)=0;

% For magnitude of normals
for i =1:N
     slope(i) = - eps* (4-8*x_c1(i));
end

angle = atand(slope);

for i = 1:N
    norm(i,1) = sind(angle(i));
    norm(i,2) = cosd(angle(i));
    n_loc(i,1) = x_c1(i);
    n_loc(i,2) = z_c1(i);
end

% Influence points
for i=1:N
    for j=1:N
        [u(i,j),w(i,j)] = vor(1,x_c1(i),z_c1(i),x_01(j),z_01(j));
        a(i,j) = dot([u(i,j) w(i,j)],[norm(i,1) norm(i,2)]);
    end
end

R = inv(a);
gam = -R*(norm(:,1));
plot (x,z);
hold on;
plot (x_01,gam/(0.05));

% Calculates the influence of vortex at (x1,z1)
function [u,w] = vor(gam,x,z,x1,z1)
    u = 0;
    w = 0;
    rx = x - x1;
    rz = z - z1;
    r = (rx^2 + rz^2);   
    v = gam/(2*pi*r);
    u = v*rz;
    w = v*(-rx);
end