close all ;
clear all ;
m = 30.51 ; % W/g
W = 299 ; B = 306 ;
Ix = 0.177 ; Iy = 3.45 ; Iz = 3.45 ;
Xcg = -1.96e-2 ; Ycg = 0 ; Zcg = 0 ;
% thruster constraints
maxXprop = 23.14369 ; % in Newtons
minXprop = -18.1423 ; % in Newtons

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

% figures
fig1 = figure() ;
ax1 = axes('Parent', fig1) ;
h1 = animatedline(ax1) ;
%{
fig2 = figure() ;
ax2 = axes('Parent', fig2) ;
h2 = animatedline(ax2) ;

fig3 = figure() ;
ax3 = axes('Parent', fig3) ;
h3 = animatedline(ax3) ;
%}
% for simulation
dt = 0.01 ;
time_steps = 10e4 ;

% initial conditions
u = 0 ; v = 0; r = 0; psi = 0 ; Xprop = 0 ; ud = 0.5 ; 
% Controller Params
Kp= 4800 ; Ki = 5000 ; error_sum = 0 ; 

for i =0:dt:100
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
   error = ud - u ;
   error_sum = error_sum + error*dt ;
   Xprop = error*Kp + error_sum*Ki;
   if Xprop > maxXprop
       Xprop = maxXprop ;
   elseif Xprop < minXprop 
       Xprop = minXprop ;
   end
   addpoints(h1,i,u) ;
  % addpoints(h2,i,v) ;
  % addpoints(h3,i,r) ;
end
