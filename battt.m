
Rthrus =0.18; % in ohm
L = 0.077;
Kthrus =0.040;

thruster_data = xlsread('T100_T_P_C.xlsx') ;
% lbf to newton conversion factor
lbf_to_N = 4.44822 ;
thrust_given = thruster_data(:,1).*lbf_to_N ;   % in Newtons
rpm_given = (2*pi/60).*thruster_data(:,4) ;     % in rps
current_given =thruster_data(:,3) ; % in Amp

[~, ind] = unique(thrust_given, 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(thrust_given, 1), ind);
% duplicate values
thrust_given(duplicate_ind) = [] ;
rpm_given(duplicate_ind) = [] ;
current_given(duplicate_ind) = [] ;
% battery discharging in V
fig1 = figure() ;
ax1 = axes('Parent', fig1) ;
h1 = animatedline(ax1) ;
E0=16;
Rbatt=0.001;
Kbatt=0.00076 ;
A1= 0.9422 ;
Q= 24*3600 ;
m = 30.51 ; % W/g
W = 299 ; B = 306 ;
Ix = 0.177 ; Iy = 3.45 ; Iz = 3.45 ;
Xcg = -1.96e-2 ; Ycg = 0 ; Zcg = 0 ;
% for x 
Xu_dot = -9.30e-001 ;
Xuu = -1.62e+000 ;
Xvr = +3.55e+001 ;
Xrr = -1.93e+000 ;

% for y  
Yv_dot = -3.55e+001 ;
Yr_dot = +1.93e+000 ;
Yvv = -1.31e+003 ;
Yrr = +6.32e-001 ;
Yur = +5.22e+000 ;
Yuv = -2.86e+001 ;

%for N 
Nr_dot = -4.88e+000 ;
Nv_dot = +1.93e+000 ;
Nvv = -3.18e+000 ;
Nrr = -9.40e+001 ;
Nur = -2.00e+000 ;
Nuv = -2.40e+001 ;
%figure
fig2 = figure() ;
ax2 = axes('Parent', fig2) ;
h2 = animatedline(ax2) ;

dt=0.001;
% initial conditions for thrust
u = 0 ; v = 0; r = 0; psi = 0 ; Xprop = 0 ; ud = 0.63 ; K= 800 ; Ki = 1000;
Kt=0.0084; D=0.076; rho=1025; omega=0;
%Initial Condition for battery
i0=0; it=0; v0=E0; w=3000*(2*pi/60) ; t = 0;
P=0; deltai=0;
depth=0;
energy=0;
e = 0 ;
vel = 0 ;
depth_vector = 0 ;
thrust = 0 ;
flag = 1 ;
prevu = 0 ;
for i =0:dt:170
    if depth == 100
        break
    end
   A = [ m-Xu_dot ,0,0 ; 0, m-Yv_dot , m*Xcg-Yr_dot ; 0, m*Xcg-Nv_dot, Iz-Nr_dot]; 
   X = [ 
       (W-B)*cos(psi) + Xuu*u*u + (Xvr+m)*v*r + (Xrr+m*Xcg)*r*r + Xprop 
       -(W-B)*sin(psi) + Yvv*v*v + Yrr*r*r + (Yur-m)*u*r + Yuv*u*v 
       -Xcg*W*sin(psi) + Nvv*v*v + Nrr*r*r + (Nur-m*Xcg)*u*r + Nuv*u*v 
       ] ;
   mat = A\X ;
   u = u + mat(1)*dt ;
   v = v + mat(2)*dt ;
   r = r + mat(3)*dt ;  
   psi = psi + r*dt ;
   e = ud - u ;
   e = e+e*dt ;
   Xprop = (ud-u)*K + Ki*e;
   if Xprop > max(thrust_given)
       Xprop = max(thrust_given) ;
   end
   %omega=(2*pi*(sqrt(Xprop/(rho*(D^4)*Kt))))/60;
   %omega = spline(thrust_given, rpm_given, Xprop) ;
   % deltai = (v0-i0.*Rthrus-Kthrus.*omega)./L  ;
   i0 = spline(thrust_given, current_given, Xprop) ;
    v0= E0-Rbatt*i0 - Kbatt*(Q/(Q-it))*(it) + A1* exp(-B*it);
   % i0=deltai*dt + i0;
    P=v0*i0;
    it=it+i0*dt;
    addpoints(h1,depth,energy) ;
    addpoints(h2,i,v0);
    depth=depth+u*dt;
    energy=energy+P*dt;
    t = t + dt ;
  % addpoints(h2,i,v) ;
  % addpoints(h3,i,r) ;
  %if prevu ~= u
  vel(flag) = u ;
  thrust(flag) = Xprop ;
  depth_vector(flag) = depth ;
  flag = flag+1 ;
  %end
  %prevu = u ;
end


