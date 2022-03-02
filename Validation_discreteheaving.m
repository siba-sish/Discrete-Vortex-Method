%Code for validation of force coefficients of a flat plate using discrete
%heaving method
 clc;
 clear;

%initialisation
n = input('input the number of time steps-');
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
normal = zeros(n,1);

%initialisation for panel discretization
p = 20;
num = p;
dx = c/p;
a1 = zeros(p+1);
rhs = zeros(p+1,1);
q = zeros(p,1);
dpre = zeros(p,1); %pressure difference between edges of panels
gammat = zeros(n,p);

%initial vortex influence
m = 0;
u_wake = 0;
w_wake = 0;
k = 0;
rhs2 = 0;
u_p = 0;
w_p = 0;
gamat1 = 0;

%Argot of program
alpha = 0;
sn = sind(alpha);
cs = cosd(alpha);
h = 0.019*c;
omet = 4.26*2*pi;
ut = 1.56*c;
dt = 0.009*c/ut;
t = 0:dt:(n-1)*dt;
norm = [sn cs];
z = -h*sin(omet*t);
wbody = h*omet*cos(omet*t);
q = sqrt(ut^2+wbody.^2);
dxw = 0.3*q*dt;
sx = -ut*t;
vortic(:,1) = (c+dxw).*cs+sx;
vortic(:,2) = -(c+dxw).*sn+z;

i = 1:p;
x_01 = (i-0.75)*dx;
x_c1 = (i-0.25)*dx;
axis tight manual
ax = gca;
v = VideoWriter('siba.free');
open(v);

%Program start
for it = 1:n
   it1 = it-1;
   cutoff = 0.0001;
   wwake = zeros(1,p);   
   if (it==1)

   else
  %    calculate wake influence for zero normal boundary condition
  %    Inner loop scans all vortices for a single collocation point at time step
  %    Outer loop finds out the total velocity at the collocation point
         it1 = it - 1;
         u1 = zeros(it1,p);
         w1 = u1;
         rx1 = u1;
         rz1 = u1;
         r = u1;
         %to find downwash and total circulation due to n-1 vortexes
         for j = 1:it1
             rx1(j,:)= x_c1+sx(it)-vortic(j,1);
             rz1(j,:)= z(it)-vortic(j,2);
             r(j,:) = rx1(j,:).^2 + rz1(j,:).^2 ;
             u1(j,:) = rz1(j,:)/(2*pi).*vortic(j,3)./r(j,:).*(1 - exp(-r(j,:)/cutoff^2));
             w1(j,:) = -rx1(j,:)/(2*pi).*vortic(j,3)./r(j,:).*(1 - exp(-r(j,:)/cutoff^2));
         end
         u = sum(u1);
         w = sum(w1);
         % downwash due to n-1 vortexes
         wwake(1,:) = u.*norm(1) + w.*norm(2);
         % total n-1 wake circulation
        rhs2 = -sum(vortic(1:it1,3));
   end
   %Calculation of wake influence coefficients
   for i=1:p
        for j=1:p   %Bound on bound
            [u(i,j),w(i,j)] = vor(x_c1(i),z(it),x_01(j),z(it),1);
            a1(i,j) = u(i,j)*norm(1) + w(i,j)*norm(2);        
        end
        [u(i,p+1),w(i,p+1)] = vor(x_c1(i),z(it),c+dxw(it)*cs,vortic(it,2),1); %vortex influence
        a1(i,p+1) = u(i,p+1)*norm(1) + w(i,p+1)*norm(2);
        rhs(i) = -(ut*norm(1)+wbody(it)*norm(2))-wwake(1,i);
   end
   a1(p+1,:) = 1;
   rhs(p+1) = rhs2; % this correct as the third row adds up.
   R = inv(a1);
   gam = R*rhs;
   vortic(it,3) = gam(p+1);
   s = sum(gam(1:p));
   %stores circulation for each time step
   for h =1:p
       gammat(it,h) = gam(h);
   end
   %  wake rollup
   if (it<1)
   else
             rx2 = zeros(it,p);
             rz2 = rx2;
             r2 = rx2;
             u4 = rx2;
             w4 = rx2;
             %influence of bound vortexes on wake vortexes
             for m = 1:it
             %influence of all bound vortices shed till n = it steps
             rx2(m,:) = vortic(m,1)-sx(it)-x_01;
             rz2(m,:) = vortic(m,2)-z(it);
             r2(m,:) = rx2(m,:).^2+rz2(m,:).^2;
             u4(m,:) = (1/2*pi).*rz2(m,:).*gammat(it,:)./r2(m,:).*(1 - exp(-r2(m,:)/cutoff^2));
             w4(m,:) = -(1/2*pi).*rx2(m,:).*gammat(it,:)./r2(m,:).*(1 - exp(-r2(m,:)/cutoff^2));
             u4(isnan(u4)) = 0;
             w4(isnan(w4)) = 0; 
             u5 = sum(u4,2);
             w5 = sum(w4,2);
             end
             
        %Determine the influence of the wake elements on each other:
        for e = 1:it
            %distance b\w all elements and current element:
            rx = vortic(e,1) - vortic(1:it,1);
            rz = vortic(e,2) - vortic(1:it,2);
            r2 = rx.^2 + rz.^2 + 0*cutoff^2;
            
            u11 = rz/(2*pi).*vortic(1:it,3)./r2.*(1 - exp(-r2/cutoff^2));
            w11 = -rx/(2*pi).*vortic(1:it,3)./r2.*(1 - exp(-r2/cutoff^2));
            
           
            u11(isnan(u11)) = 0;
            w11(isnan(w11)) = 0;
           
            u6 = sum(u11);
            w6 = sum(w11);
            u5(e) = u5(e) + u6;%final x velocities of vortices
            w5(e) = w5(e) + w6;%final w velocities of vortices
        end
       vortic(1:it,1) = vortic(1:it,1) + u5(1:it)*dt;%final positions
       vortic(1:it,2) = vortic(1:it,2) + w5(1:it)*dt;
   end     
    %Aerodynamic loads
    gamat1 = 0;
    if (it==1)
        gamat1 = 0;
    end
        q = 0.5*rho*ut*ut;
        gammanet = sum(gammat(it,:));
        dgamdt = (gammanet-gamat1)/dt;
        if (it==1)
           dgamdt1 = (gammat(it,:) - zeros(1,p))/dt;
        else
        dgamdt1 = (gammat(it,:) - gammat(it-1,:))/dt;
        end
        gamat1 = gammanet;    
        for i = 1:it
            for m = 1:p
                [uvel_1(i,m),wvel_1(i,m)] = vor(x_c1(m)+sx(it),z(it),vortic(i,1),vortic(i,2),vortic(i,3));
            end
        end
        u_v = sum(uvel_1);
        w_v = sum(wvel_1);
        wwake1 = u_v*norm(1) + w_v*norm(2);
        %Drag and Pressure difference calculation
         for i = 1:p
             drag(i) = rho*(wwake1(i)*gam(i)+(sum(dgamdt1(1,1:i)*dx*norm(1))));
             dpre(i) = rho*(((ut+u_v(i))*norm(2)+((wbody(it)+w_v(i))*-norm(1)))*gam(i)/dx + sum(dgamdt1(1,1:i)));
             lift(i) = dpre(i)*dx*norm(2);
         end
        cl(it) = sum(lift)*c/q;
        cd(it) = sum(drag)*c/q;
    normal(it) = cl(it)*cs + cd(it)*sn;
    plot(vortic(1:it,1)/c,vortic(1:it,2)/c,'.');
    axis([-5 2 -1.5 1.5])
    grid on;
    xlabel('x');
    ylabel('z');
    frame = getframe(gcf);
    writeVideo(v,frame);
end
 close(v)
 time = ut*t/c;
 figure
 plot(vortic(:,1)/c,vortic(:,2)/c,'.')
 hold on;
 plot(sx,z,'-r')
 xlabel('x')
 ylabel('y')
 figure
 plot(time(10:n)/c,cl(10:n));
 xlabel('t')
 ylabel('Lift Coefficient')
 figure
 plot(time(10:n)/c,cd(10:n));
 xlabel('t')
 ylabel('Drag Coefficient')
 figure
 plot(-sx,z,'.');
 hold on;
 plot((time(10:n))/c,normal(10:n));
 xlabel('t')
 ylabel('cn')

%calculates the influence of vortex at (x1,z1)
function [u,w] = vor(x,z,x1,z1,gam)
    cutoff = 0.0001;
    u = 0;
    w = 0;
    rx = x - x1;
    rz = z - z1;
    r = rx^2 + rz^2;   
    v = gam/(2*pi*r);
    u = v*rz*(1 - exp(-r/cutoff^2));
    w = v*(-rx)*(1 - exp(-r/cutoff^2));
end