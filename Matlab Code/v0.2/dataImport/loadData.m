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

%Define Buses
run('loadAnalysisCasesBus.m')
run('loadPlantDataBus.m')
run('loadEnvDataBus.m')


%% CHECK DIMENSIONS
run('dataValidations.m')

% busInfo = Simulink.Bus.createObject('untitled1', 'untitled1/analysisCasesBus','object')

