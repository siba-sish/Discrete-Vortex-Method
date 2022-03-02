%Transient lumped vortex method for sudden accleration of a flat plate
clc;
clear;

%initialisation
n = 40;
rho = 1.225;
ut = 50;
c = 1;
uk = 0;
wk = 0;
it = 0;
m = 0;
x = linspace(0,1,n);
vortic = zeros(n,3);
gammat = 0;
dgamat = 0;
gamat1 = 0;
drag = zeros(n,1);
time = zeros(n,1);
ratio = zeros(n,1);
uw = zeros(n,2);

%Argot of program
alpha = 5;
sn = sind(alpha);
cs = cosd(alpha);
dt = c*0.25 / ut;
t = -dt;
dxw = 0.3*ut*dt;

% plot(x,-tand(alpha)*x);
% hold on;

%Program start The inclination accounted for in collocation and vortex
%points

for it = 1:n
   t = t+dt;
   %path of origin
   sx = -ut*t;
   sz = 0;
%    plot(x+sx,-tand(alpha)*x);
%    hold on;
   %shedding of wake points
   vortic(it,1) = (c+dxw)*cs+sx;
   vortic(it,2) = -(c+dxw)*sn+sz;
%    plot(vortic(it,1),vortic(it,2),'o');
%    hold on;
   %calculate momentary vortex strength of wing and wake vortices
   a = -1/(pi*c);
   b = 1/(2*pi*(c/4+dxw));
   rhs2 = 0;

   if (it==1)
       continue
   else
       it1 = it-1;
       %calculate wake influence for zero normal boundary condition
       xx1 = 0.75*c*cs+sx;
       zz1 = -0.75*c*sn+sz;
%        plot(xx1,zz1,'.');
%        hold on;
       u = 0;
       w = 0;
       for i = 1:it1
         [u1,w1] = vor(xx1,zz1,vortic(i,1),vortic(i,2),vortic(i,3));        
         u = u + u1;
         w = w + w1;
       end
       wwake = u*sn+w*cs;
       %calculation of rhs
       for i = 1:it1
           rhs2 = rhs2-vortic(i,3);
       end
    end   
    rhs1 = -ut*sn-wwake;
    %solution based on algebraic solution of two equations for gmmat and the latest wake vortex strength vortic(it,3))
    vortic(it,3) = 1/((b/a)-1)*((rhs1/a)-rhs2);
    gammat = rhs2-vortic(it,3)
 
    %wake rollup
    if (it<1)
        continue
    else
        for i = 1:it
            xx1 = 0.25*c*cs+sx;
            zz1 = -0.25*c*sn+sz;
            [u,w] = vor(vortic(i,1),vortic(i,2),xx1,zz1,gammat);
            uk = 0;
            wk = 0;
            for m = 1:it1
                if (m==i)
                    continue;
                else
                    [u1,w1] = vor(vortic(i,1),vortic(i,2),vortic(m,1),vortic(m,2),vortic(m,3));        
                    uk = uk + u1;
                    wk = wk + w1;  
                end
            end
            u = u+uk;
            w = w+wk;
            uw(i,1) = vortic(i,1)+u*dt;
            uw(i,2) = vortic(i,2)+w*dt;
            vortic(i,1) = uw(i,1);
            vortic(i,2) = uw(i,2);         
        end
    end

    %Aerodynamic loads
    if (it==1)
        gamat1=0;
        continue
    else
        q = 0.5*rho*ut*ut;
        dgamdt = (gammat-gamat1)/dt;
        gamat1 = gammat;
        %Calculate wake induced downwash
        xx1 = 0.75*c*cs+sx;
        zz1 = -0.75*c*sn+sz;
        u = 0;
        w = 0;
        for i = 1:it1
           [u1,w1] = vor(xx1,zz1,vortic(i,1),vortic(i,2),vortic(i,3));        
           u = u + u1;
           w = w + w1;      
        end        
        ww = u*sn+w*cs;
        l = rho*(ut*gammat+dgamdt*c);
        l1 = rho*(ut*gammat);
        d = rho*(-ww*gammat+dgamdt*c*sn);
        cl = l*c/q;
        cd = d*c/q;
        %output
        clt = cl/(2*pi*sn);
        gam1 = gammat/(pi*ut*c*sn);
        sx1 = sx-ut*dt;
        ratio(i) = l/l1;
%         disp(cd)
        drag(i) = cd;
        time(i) = ut*t/c;
    end
end
%  drag;
%  vortic(1:n,:)
%  gammat;
%  sum(vortic(:,3));
 plot(time(2:n-1),drag(2:n-1),'--');
 hold on;
% plot(vortic(:,1),vortic(:,2));

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