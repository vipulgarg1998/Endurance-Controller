close all ;
clear ;
%%%%%%%%%%%%%%%%%%%% Vehicle constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
available_energy = 50000 ;  % energy stored in battery (in joules)
dives_req = 6    ;          % number of dives to be performed 
high_speed_operation = false;% operate the vehicle at higher speeds
dives_left = dives_req ;    % dives left

%%%%%%%%%%%%%%%%%%%% get the data for velocity and energy %%%%%%%%%%%%%%%%%
[min_velocity, min_energy, index, velocity_of_dives, energy_of_dives, load] = compute_power() ;
energy_of_dives = energy_of_dives + load;

%%%%%%%%%%%%%%%%%%%% Simulate REMUS                     
net = Remus_model(min_velocity) ;
[min_velocity, min_energy, index, velocity_of_dives, energy_of_dives, load] = compute_power_NN(net) ;

%%%%%%%%%%%%%%%%%%%% Initialisation %%%%%%%%%%%%%%%%%%%%%
vel = zeros(dives_req,1) ;
ener = zeros(dives_req,1) ;
for dive = 1:1:dives_req
    expected_energy_consumption = max(available_energy/dives_left,min_energy) ;
    if high_speed_operation
        vel(dive) = spline(energy_of_dives(index:end),velocity_of_dives(index:end),expected_energy_consumption);
    else
        vel(dive) = spline(energy_of_dives(1:index),velocity_of_dives(1:index),expected_energy_consumption);
    end
    energy_consumed = 0.9*expected_energy_consumption + 0.2*expected_energy_consumption*rand();
    available_energy = available_energy - energy_consumed ;
    ener(dive) = available_energy ;
    dives_left = dives_left - 1;
    energy_of_dives = energy_of_dives + energy_consumed - expected_energy_consumption ;
end