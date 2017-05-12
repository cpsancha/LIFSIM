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

%Input of landing gear properties
run('loadLandingGearData.m')

%Input of FCS Data
run('loadFCSData.m')

%Input of Sensors and Actuators properties
run('loadSensorsActuators.m')

%Define Buses
run('loadAnalysisCasesBus.m')
run('loadPlantDataBus.m')
run('loadEnvDataBus.m')
run('loadSensorsBus.m')
run('loadCommandBus.m')
run('loadThrottleBus.m')
run('loadActuatorsBus.m')


%% CHECK DIMENSIONS
run('dataValidations.m')


