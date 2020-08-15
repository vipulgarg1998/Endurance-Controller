# Endurance-Controller
A custom package to determine the optimum velocity and control the endurance of the Autonomous Vertical Profiler.

# MATLAB code
change the MATLAB directory to Endurance-Controller and run the following commands
```
addpath('mat_files\');
main
```
# Output
plot velocity v/s energy
``` 
plot(velocity_of_dives(:), energy_of_dives(:)) ;
```
display optimum velocity and corresponding energy consumed during a dive
```
display(min_velocity) ; 
display(min_energy)   ;
```
display the vel0city and available energy for each dive 
```
display(vel) ;
display(ener);
```
