% Transient discrete vortex method for sudden accleration of a flat plate
% Flat plate in a specific angle of attack and velocity lateral
clc;
clear;

%initialisation
n = 40;
rho = 1.225;
c = 1;
it = 0;
vortic = zeros(n,3);
gam = zeros(n,1);
drag = zeros(n,1);
lift = zeros(n,1);
time = zeros(n,1);
uvel_1 = zeros(n,1);
wvel_1 = zeros(n,1);
uvel = zeros(n,1);
wvel = zeros(n,1);
cd = zeros(n,1);
cl = zeros(n,1);
u_w = zeros(n,1);
w_w = zeros(n,1);
u_v = zeros(n,1);
w_v = zeros(n,1);
ubw = zeros(n,1);
wbw = zeros(n,1);
uww = zeros(n,1);
www = zeros(n,1);
wwake1 = zeros(n,1);

%initialisation for panel discretization
p = 100;
dx = c/p;
x_01 = zeros(p,1);
x_c1 = zeros(p,1);
z_01 = zeros(p,1);
z_c1 = zeros(p,1);
z_01re = zeros(p,1);
z_c1re = zeros(p,1);
a1 = zeros(p+1);
rhs = zeros(p+1,1);
q = zeros(p,1);
dpre = zeros(p,1); %pressure difference between edges of panels
dp = c/p;
gammat = zeros(n,p);

%initial vortex influence
m = 0;
u_wake = 0;
w_wake = 0;
k = 0;
s = 0;
u_p = 0;    
w_p = 0;
gamat1 = 0;

%Argot of program
alpha = 5;
sn = sind(alpha);
cs = cosd(alpha);
ut = 50;
zt = 0;
dt = c*0.25/ut;
t = -dt;
dxw = 0.3*ut*dt;
norm = [sn cs];

i = 1; 
x_01(i) = 0.25*dx ;
x_c1(i) = 0.75*dx ;

z_01(i) = -tand(alpha)*x_01(i) ;
z_c1(i) = -tand(alpha)*x_c1(i) ;
z_01re(i) = z_01(i);
z_c1re(i) = z_c1(i);

for  i =2:p
   x_01(i) = dx + x_01(i-1);
   x_c1(i) = dx + x_c1(i-1);
   z_01(i) = -tand(alpha)*x_01(i) ;
   z_c1(i) = -tand(alpha)*x_c1(i) ;
   z_01re(i) = z_01(i);
   z_c1re(i) = z_c1(i);
end
%    plot(x_01,z_01);
%    hold on;

%Calculation of bound infulence coefficients
for i=1:p
     for j=1:p   %Bound on bound
        [u(i,j),w(i,j)] = vor(x_c1(i),z_c1(i),x_01(j),z_01(j),1);
        a1(i,j) = u(i,j)*norm(1) + w(i,j)*norm(2);
     end 
end

%Program start
for it = 1:n
   t = t+dt;
   %path of origin
   sx = -ut*t;
   sz = 0;
   vortic(it,1) = (c+dxw)*cs+sx;
   vortic(it,2) = -(c+dxw)*sn+sz;
   i = 1;
   x_01(i) = 0.25*dx + sx;
   x_c1(i) = 0.75*dx + sx;
  
   % Coordinates of required points in panels
   if (it==1)
        continue;
   else
        z_01(1) = z_01re(1);
        z_c1(1) = z_c1re(1);
        for i = 2:p
            x_01(i) = dx + x_01(i-1);
            x_c1(i) = dx + x_c1(i-1);
            z_01(i) = 0 + z_01re(i);
            z_c1(i) = 0 + z_c1re(i);
        end
   end
   if (it==1)
       continue
   else
       it1 = it-1;
  %    calculate wake influence for zero normal boundary condition
  %    Inner loop scans all vortices for a single collocation point at time step
  %    Outer loop finds out the total velocity at the collocation point
        for i = 1:p
            for m = 1:it1 %it1 due to unknown latest wake vortex
                [uvel_1(m),wvel_1(m)] = vor(x_c1(i),z_c1(i),vortic(m,1),vortic(m,2),vortic(m,3));
            end
            u_w(i) = sum(uvel_1);
            w_w(i) = sum(wvel_1);
        end
        wwake = sum(u_w)*sn + sum(w_w)*cs;
        rhs2 = -sum(vortic(1:it1,3));
   end
   %Supposedly correct downwash
%    u_w
%    w_w
   %Calculation of wake influence coefficients
   for i=1:p       %Wake on bound
        [u(i,p+1),w(i,p+1)] = vor(x_c1(i),z_c1(i),vortic(it,1),vortic(it,2),1); %vortex influence
        a1(i,p+1) = u(i,p+1)*norm(1) + w(i,p+1)*norm(2);
        rhs(i) = -((ut+u_w(i))*norm(1)+(zt+w_w(i))*norm(2));
   end
   a1(p+1,:) = 1;
   rhs(p+1) = rhs2; % this correct as the third row adds up.
   R = inv(a1);
   gam = R*rhs;;
   vortic(it,3) = gam(p+1);
   s = sum(gam(1:p));
   %stores circulation for each time step
   for h =1:p
       gammat(it,h) = gam(h);
   end

    %Aerodynamic loads
    if (it==1)
        continue;
    else
        q = -0.5*rho*ut*ut;
        gammanet = sum(gammat(it,:));
        dgamdt = (gammanet-gamat1)/dt;
        if (it==1)
           dgamdt1 = (gammat(it,:) - zeros(1,num))/dt;
        else
        dgamdt1 = (gammat(it,:) - gammat(it-1,:))/dt;
        end
        gamat1 = gammanet;
        circulation(it) = gammanet;
        
        %Calculate wake induced downwash    
        for i = 1:p
            for m = 1:it1
                [uvel_1(m),wvel_1(m)] = vor(x_c1(i),z_c1(i),vortic(m,1),vortic(m,2),vortic(m,3));
            end
                u_v(i) = sum(uvel_1);
                w_v(i) = sum(wvel_1);
                wwake1(i) = u_v(i)*norm(1)+w_v(i)*norm(2);
        end
        ww = sum(wwake1);
        %Drag and Pressure difference calculation
         for i = 1:p
             drag(i) = rho*(wwake1(i)*gam(i)-(sum(dgamdt1(1,1:i)*dp*norm(1))));
             dpre(i) = rho*(((ut+u_v(i))*norm(2)+((zt+w_v(i))*-norm(1))*gam(i)/dp) + sum(dgamdt1(1,1:i)));
             lift(i) = dpre(i)*dp*norm(2);
         end
        cl(it) = sum(lift)*c/q;
        cd(it) = sum(drag)*c/q;
        time(it) = ut*t/c;
    end
    
    %wake rollup
     for k = 1:it
         for i = 1:p    %Bound on wake        
             [ubw(i),wbw(i)] = vor(vortic(k,1),vortic(k,2),x_01(i),z_01(i),gam(i));
         end  
         for m = 1:k-1  %Wake on wake
             [uww(m),www(m)] = vor(vortic(k,1),vortic(k,2),vortic(m,1),vortic(m,2),vortic(m,3));
         end
         u_wake = sum(ubw)+sum(uww);
         w_wake = sum(wbw)+sum(www);
         u_p = vortic(k,1)+u_wake*dt;
         w_p = vortic(k,2)+w_wake*dt;
         vortic(k,1) = u_p;
         vortic(k,2) = w_p;
     end
end
 plot(time(2:n),cd(2:n));

%calculates the influence of vortex at (x1,z1)
function [u,w] = vor(x,z,x1,z1,gam)
    u = 0;
    w = 0;
    rx = x - x1;
    rz = z - z1;
    r = sqrt(rx^2 + rz^2);   
    v = gam/(2*pi*(r^2));
    u = v*rz;
    w = v*(-rx);
end