%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to fill all the LIBIS data into the structure LD. Each subscript 
% can be edited individually to define different parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Remove residual information
clear all %TO BE REMOVED
clear LD

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

%Define Variants
Variants.Motor1=1;
Variants.Motor2=1;
Variants.Motor3=1;
Variants.Motor4=1;
Variants.Motor5=1;

% Add enum structure for the Variants
 Simulink.defineIntEnumType('Variants',...
     {'Motor1','Motor2','Motor3','Motor4','Motor5'},[1;1;1;1;1]);

%% CHECK DIMENSIONS
run('dataValidations.m')


% busInfo = Simulink.Bus.createObject('untitled1', 'untitled1/analysisCasesBus','object')

