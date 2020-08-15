function [min_velocity, min_energy, index, velocity_of_dives, energy_of_dives, load] = compute_power_NN(net)

%%%%%%%%%%%%%%%%% read data from CTD dataset %%%%%%%%%%%%%%%%%%%%%
[depth, density_data] = load_CTD_data('NPEO2000_CTD1', 'mat_files\seawater\') ;

%%%%%%%%%%%%%% read data from Thruster Dataset %%%%%%%%%%%%%%%%
[thrust_given, power_given] = load_thruster_data('mat_files\T100_T_P_C.xlsx') ;

%%%%%%%%%%%%%%%%%%%%% system contants %%%%%%%%%%%%%%%%%%%%%%%%%%
body_mass = 30.51;% 30.51 %18.25; % 50
body_volume = 0.0305;% 0.0305 %0.018;    % 0.049
drag_coefficient = 0.27; % 0.2
reference_area = 0.028; % 0.1
g = 9.8 ;
G_F = body_mass*g ;
Scientific_Sensors =  4.26 ;
Vehicle_Sensors = 3.382 ;
Electronics = 2.52 ;

%%%%%%%%%%%%%%%%%%%% Simulation Parameters %%%%%%%%%%%%%%%%
desired_depth = 100;
no_of_steps = 100;
initinal_velocity = 0.1 ;
increment_in_velocity = 0.01 ;
final_velocity = 2.4 ;
velocity_of_dives = (initinal_velocity:increment_in_velocity:final_velocity) ;
no_of_observations = size(velocity_of_dives,1) ;
energy_of_dives = zeros(no_of_observations,1) ;
dive_number = 1 ;
depth_step = 1 ;
%%%%%%%%%%%%%%%%%%%%%% main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
energy_record = zeros(no_of_steps,1) ;
thrust_record = zeros(no_of_steps,1) ;
for vel = velocity_of_dives
    energy_of_dives(dive_number) = 0 ;
    time_step = (desired_depth/no_of_steps)/vel ;
    for dep = linspace(desired_depth/no_of_steps,desired_depth,no_of_steps)
        rho = spline(depth,density_data,dep) ;
        B_F = rho*body_volume*g ;
        if vel < -1
            D_F = drag_coefficient*reference_area*rho*(vel)/2 ;
        else
            D_F = -net(vel) ;
        end
        T_F = B_F + D_F - G_F ;
        power_consumed = spline(thrust_given, power_given, T_F) ;
        energy_of_dives(dive_number) = energy_of_dives(dive_number) + (power_consumed + Scientific_Sensors + Vehicle_Sensors + Electronics)*time_step ;
        thrust_record(depth_step) = T_F ;
        energy_record(depth_step) = energy_of_dives(dive_number);
    end
    if dive_number ==63 
        %plot(ax2,linspace(desired_depth/no_of_steps,desired_depth,no_of_steps),thrust_record) ;
        %plot(ax3,linspace(desired_depth/no_of_steps,desired_depth,no_of_steps),blah) ;
    end
    dive_number = dive_number + 1;
    
end
%plot(ax1,velocity_of_dives(:), energy_of_dives(:)) ;
%%%%%%%%%%%%%%%% results %%%%%%%%%%%%%%%%%%
[min_energy, index] = min(energy_of_dives) ;
min_velocity = min(velocity_of_dives(index)) ;
result = sprintf('the minimum energy is %f Joules for %f m/s',min_energy, min_velocity) ;
display(result) ;
%%%%%%%%%%%%%%%% upward motion %%%%%%%%%%%%
depth_step = 1 ;
u = 0 ;
time = 0 ;
for dep = desired_depth:-depth_step:1
    rho = spline(depth,density_data,dep) ;
    B_F = rho*body_volume*g ;
    D_F = drag_coefficient*reference_area*rho*(u^2)/2 ;
    a = (B_F-D_F-G_F)/(body_mass) ;
    v = sqrt(2*a*depth_step + u^2) ;
    time = time + (v - u)/a ;
    u = v ;
end
load = (Scientific_Sensors + Vehicle_Sensors + Electronics)*time ;