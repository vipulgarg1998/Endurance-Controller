fig1 = figure() ;
ax1 = axes('Parent', fig1) ;
axis(ax1,'equal') ;
fig3 = figure() ;
ax3 = axes('Parent', fig3) ;
hold(ax3, 'on') ;

%%%%%%%%%%%%%%%%% read data from CTD dataset %%%%%%%%%%%%%%%%%%%%%
data = importdata('NPEO2000_CTD1') ;
CTD_data = data.data(:,:) ;
CTD_data_size = size(CTD_data,1) ;

density_data = zeros(CTD_data_size,1) ;
depth = CTD_data(:,1) ;
T_F = zeros(CTD_data_size,1) ;
time = zeros(CTD_data_size,1) ;
current_dir = cd('mat_files\seawater\') ;
files = dir('*.m') ;
names = cell(size(files)) ;
parfor i = 1:max(size(files))
    names{i} = files(i).name ;
end
fig1 = figure() ;
ax1 = axes('Parent', fig1) ;
axis(ax1,'equal') ;
fig3 = figure() ;
ax3 = axes('Parent', fig3) ;
hold(ax3, 'on') ;

%%%%%%%%%%%%%%%%% read data from CTD dataset %%%%%%%%%%%%%%%%%%%%%
data = importdata('NPEO2000_CTD1') ;
CTD_data = data.data(:,:) ;
CTD_data_size = size(CTD_data,1) ;

density_data = zeros(CTD_data_size,1) ;
depth = CTD_data(:,1) ;
T_F = zeros(CTD_data_size,1) ;
time = zeros(CTD_data_size,1) ;
current_dir = cd('mat_files\seawater\') ;
files = dir('*.m') ;
names = cell(size(files)) ;
parfor i = 1:max(size(files))
    names{i} = files(i).name ;
end
poolobj = gcp;
addAttachedFiles(poolobj,names') ;
cd(current_dir) ;
parfor i = 1:max(size(CTD_data))
    % density = f( salinity, temperature, pressure)
    density_data(i) = sw_dens(CTD_data(i,5), CTD_data(i,3), CTD_data(i,2)) ;
end

%%%%%%%%%%%%%% read data from Thruster Dataset %%%%%%%%%%%%%%%%
thruster_data = xlsread('mat_files\T100_T_P_C.xlsx') ;
% lbf to newton conversion factor
lbf_to_N = 4.44822 ;
thrust_given = thruster_data(:,1).*lbf_to_N ;   % in Newtons
power_given = thruster_data(:,2) ;              % in Watts

[~, ind] = unique(thrust_given, 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(thrust_given, 1), ind);
% duplicate values
thrust_given(duplicate_ind) = [] ;
power_given(duplicate_ind) = [] ;

%%%%%%%%%%%%%%%%%%%%% system contants %%%%%%%%%%%%%%%%%%%%%%%%%%
body_mass = 50;
body_volume = 0.049;
drag_coefficient = 0.2;
reference_area = 0.1;
desired_depth = 100;
no_of_steps = 100;
initinal_velocity = 0.01 ;
increment_in_velocity = 0.01 ;
final_velocity = 1.5 ;
no_of_observations = size(initinal_velocity:increment_in_velocity:final_velocity,1) ;
g = 9.8 ;
G_F = body_mass*g ;
energy_of_dives = zeros(no_of_observations,1) ;
velocity_of_dives = zeros(no_of_observations,1) ;
dive_number = 1 ;
%%%%%%%%%%%%%%%%%%%%%% main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blah = zeros(100,1) ;
for vel = initinal_velocity:increment_in_velocity:final_velocity
    energy_of_dives(dive_number) = 0 ;
    velocity_of_dives(dive_number) = vel ;
    time_step = (desired_depth/no_of_steps)/vel ;
    for dep = linspace(desired_depth/no_of_steps,desired_depth,no_of_steps)
        rho = spline(depth,density_data,dep) ;
        B_F = rho*body_volume*g ;
        D_F = drag_coefficient*reference_area*rho*(vel^2)/2 ;
        T_F = B_F + D_F - G_F ;
        power_consumed = spline(thrust_given, power_given, T_F) ;
        
        energy_of_dives(dive_number) = energy_of_dives(dive_number) + power_consumed*time_step ;

        blah(dep) = energy_of_dives(dive_number);
    end
    if dive_number ==39 
        display(size(blah)) ;
        plot(ax3,linspace(desired_depth/no_of_steps,desired_depth,no_of_steps),blah) ;
    end
    dive_number = dive_number + 1;
    
end
plot(ax1,velocity_of_dives(:), energy_of_dives(:)) ;
%%%%%%%%%%%%%%%
depth_step = 1 ;
[min_dive, index] = min(energy_of_dives) ;
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

poolobj = gcp;
addAttachedFiles(poolobj,names') ;
cd(current_dir) ;
parfor i = 1:max(size(CTD_data))
    % density = f( salinity, temperature, pressure)
    density_data(i) = sw_dens(CTD_data(i,5), CTD_data(i,3), CTD_data(i,2)) ;
end

%%%%%%%%%%%%%% read data from Thruster Dataset %%%%%%%%%%%%%%%%
thruster_data = xlsread('mat_files\T100_T_P_C.xlsx') ;
% lbf to newton conversion factor
lbf_to_N = 4.44822 ;
thrust_given = thruster_data(:,1).*lbf_to_N ;   % in Newtons
power_given = thruster_data(:,2) ;              % in Watts

[~, ind] = unique(thrust_given, 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(thrust_given, 1), ind);
% duplicate values
thrust_given(duplicate_ind) = [] ;
power_given(duplicate_ind) = [] ;

%%%%%%%%%%%%%%%%%%%%% system contants %%%%%%%%%%%%%%%%%%%%%%%%%%
body_mass = 50;
body_volume = 0.049;
drag_coefficient = 0.2;
reference_area = 0.1;
desired_depth = 100;
no_of_steps = 100;
initinal_velocity = 0.01 ;
increment_in_velocity = 0.01 ;
final_velocity = 1.5 ;
no_of_observations = size(initinal_velocity:increment_in_velocity:final_velocity,1) ;
g = 9.8 ;
G_F = body_mass*g ;
energy_of_dives = zeros(no_of_observations,1) ;
velocity_of_dives = zeros(no_of_observations,1) ;
dive_number = 1 ;
%%%%%%%%%%%%%%%%%%%%%% main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blah = zeros(100,1) ;
for vel = initinal_velocity:increment_in_velocity:final_velocity
    energy_of_dives(dive_number) = 0 ;
    velocity_of_dives(dive_number) = vel ;
    time_step = (desired_depth/no_of_steps)/vel ;
    for dep = linspace(desired_depth/no_of_steps,desired_depth,no_of_steps)
        rho = spline(depth,density_data,dep) ;
        B_F = rho*body_volume*g ;
        D_F = drag_coefficient*reference_area*rho*(vel^2)/2 ;
        T_F = B_F + D_F - G_F ;
        power_consumed = spline(thrust_given, power_given, T_F) ;
        
        energy_of_dives(dive_number) = energy_of_dives(dive_number) + power_consumed*time_step ;

        blah(dep) = energy_of_dives(dive_number);
    end
    if dive_number ==39 
        display(size(blah)) ;
        plot(ax3,linspace(desired_depth/no_of_steps,desired_depth,no_of_steps),blah) ;
    end
    dive_number = dive_number + 1;
    
end
plot(ax1,velocity_of_dives(:), energy_of_dives(:)) ;





    
    