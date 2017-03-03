%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to fill all the LIBIS data into the structure LD. Each subscript 
% can be edited individually to define different parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Remove residual information
clear all %TO BE REMOVED
clear LD
tic


%Input of aerodynamic data
run('loadStabilityData.m')

%Input of inertia properties
run('loadInertiaData.m')

%Input of propulsion properties
run('loadPropulsionData.m')

%Define Buses
run('loadAnalysisCasesBus.m')
run('loadPlantDataBus.m')
run('loadEnvDataBus.m')


%% CHECK DIMENSIONS
run('dataValidations.m')

toc
% busInfo = Simulink.Bus.createObject('untitled1', 'untitled1/analysisCasesBus','object')

