%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to fill all the LIBIS data into the structure LD. Each subscript 
% can be edited individually to define different parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Remove residual information
clear LD

%Input of aerodynamic data
run('loadStabilityData.m')

%Input of inertia properties
run('loadInertiaData.m')

%Input of propulsion properties
run('loadPropulsionData.m')



%% CHECK DIMENSIONS
run('dataValidations.m')

