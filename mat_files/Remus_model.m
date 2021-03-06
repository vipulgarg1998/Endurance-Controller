function net = Remus_model(ud)
m = 30.51 ; % W/g
W = 299 ; B = 306 ;
Ix = 0.177 ; Iy = 3.45 ; Iz = 3.45 ;
Xcg = -1.96e-2 ; Ycg = 0 ; Zcg = 0 ;
% thruster constraints
maxXprop = 23.14369 ; % in Newtons
minXprop = -18.1423 ; % in Newtons

% for x 
Xu = -13.5 ;
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

fig2 = figure() ;
ax2 = axes('Parent', fig2) ;
h2 = animatedline(ax2) ;
%{
fig3 = figure() ;
ax3 = axes('Parent', fig3) ;
h3 = animatedline(ax3) ;
%}
% for simulation
dt = 0.001 ;
flag = 1 ;
inputs = [] ;
outputs = [] ;
% initial conditions
u = 0 ; v = 0; r = 0; psi = 0 ; Xprop = 0 ; % ud = 0.5 ; 
% Controller Params
Kp = 96.6; Ki = 29.4; Kd = 4.24 ; prev_error = 0 ; error_sum = 0 ; 

%{
Xu = -13.5 ; 
sys = zpk([], -Xu/(m-Xu_dot), 1/(m-Xu_dot)) ;
[C_pi,info] = pidtune(sys,'PID') ;
T_pi = feedback(C_pi*sys, 2);
step(T_pi)
%}
for i =0:dt:100
   % A*mat = X
   A = [ m-Xu_dot ,0,0 ; 0, m-Yv_dot , m*Xcg-Yr_dot ; 0, m*Xcg-Nv_dot, Iz-Nr_dot]; 
   X = [ 
       (W-B)*cos(psi) + Xu*u + (Xvr+m)*v*r + (Xrr+m*Xcg)*r*r + Xprop 
       -(W-B)*sin(psi) + Yvv*v*v + Yrr*r*r + (Yur-m)*u*r + Yuv*u*v 
       -Xcg*W*sin(psi) + Nvv*v*v + Nrr*r*r + (Nur-m*Xcg)*u*r + Nuv*u*v 
       ] ;
   mat = A\X ;
   drag = Xu*u ;
   force = mat(1) ;
   u = u + mat(1)*dt ;%+ 0.001*rand();
   v = v + mat(2)*dt ;
   r = r + mat(3)*dt ;
   psi = psi + r*dt ;
   error = ud - u ;
   error_sum = error_sum + error ;
   error_dif = error - prev_error ;
   prev_error = error ;
   Xprop = error*Kp + error_sum*dt*Ki + (Kd/dt)*error_dif;
   if Xprop > maxXprop
       Xprop = maxXprop ;
   elseif Xprop < minXprop 
       Xprop = minXprop ;
   end
   inputs(flag,:) = [u, drag, W-B] ;
   outputs(flag,:) = [X(1), force] ;
   addpoints(h1,i, u) ; 
   addpoints(h2,i,Xprop) ;
  % addpoints(h3,i,r) ;
  flag = flag + 1 ;
end
net = newfit(inputs', outputs', [5 5 5]);
net.performFcn = 'mae';
net.trainFcn = 'trainlm' ;
net = train(net, inputs', outputs') ;
output = net(inputs') ;
perf = perform(net,output,outputs') ;
%{
[~, ind] = unique(u_data, 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(u_data, 1), ind);
% remove duplicate values
u_data(duplicate_ind) = [] ;
x_data(duplicate_ind) = [] ;
t_data(duplicate_ind) = [] ;
%plot(t_data,u_data) ;
%}
end
