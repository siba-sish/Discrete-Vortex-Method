% Code for flow physics analysis of a thin flat airfoil using the Discrete Vortex Method
clc;
clear;

%Variable declaration and initialisation
alpha = 10;
q = 1;
N = 20;
c=1;
eps = 0.1 *c;
rho = 1.225;
x = linspace(0,c,N+1);
dx = x(3)-x(2);
x_01 = zeros(N,1);
x_c1 = zeros(N,1);
z_01 = zeros(N,1);
z_c1 = zeros(N,1);
a = zeros(N);
V = zeros(N,1);
norm1 = zeros(N,1);
norm2 = zeros(N,1);
u_inf = q*cosd(alpha);
w_inf = q*sind(alpha);

syms f(x)
y = f(x);
y = 4*eps*(x/c)*(1-(x/c));
df = @(x) (4/c)*eps*(1-2*(x/c));

x_01(1) = 0.25*dx;
x_c1(1) = 0.75*dx;

% Coordinates of required points in panels
for  i =2:size(x_01,1);
   x_01(i) = dx + x_01(i-1);
   x_c1(i) = dx + x_c1(i-1);
   z_01(i) = 4*x_01(i)*eps*(1-x_01(i));
   z_c1(i) = 4*x_c1(i)*eps*(1-x_c1(i));
end
% For magnitude of normals
norm1 = -df(x_c1)./sqrt(1+(df(x_c1).^2));
norm2 = 1./sqrt(1+(df(x_c1).^2));

% Influence points
for i=1:N
    for j=1:N
        [u(i,j),w(i,j)] = vor(1,x_c1(i),z_c1(i),x_01(j),z_01(j));
        a(i,j) = u(i,j)*norm1(i,1) + w(i,j)*norm2(i,1);
    end
    V(i) = -(u_inf*norm1(i,1)+w_inf*norm2(i,1));
end

R = inv(a);
gam =  R*V;
L = rho*q*gam;
cL = L/(0.5 * rho * q^2);
disp(sum(cL))
scatter(x_01,2*gam*(N/10)/(eps*q));
hold on;
x1 = linspace(0.01,1,N);

for i=1:N
    del_cp(i) = 4 * (((c-x1(i))/x1(i))^0.5)*sind(alpha) +32 *(eps/c)*(x1(i)*(1-x1(i)))^0.5;
end
plot (x1,del_cp);

%Calculates the influence of vortex at (x1,z1)
function [u,w] = vor(gam,x,z,x1,z1)
    u = 0;
    w = 0;
    rx = x - x1;
    rz = z - z1;
    r = sqrt(rx^2 + rz^2);   
    v = gam/(2*pi*(r^2));
    u = v*rz;
    w = v*(-rx);
end