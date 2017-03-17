%% startVars.m - Initialize variables
% This script initializes variables and buses required for the model to
% work. Mask block parameters are defined by structures that define the
% location of the block, ie. If the Initial Condition parameter is located
% under Vehicle/Nonlinear/Integrator the variable is set to
% Vehicle.Nonlinear.Integrator.initialCondition = 0;

%   Copyright 2013-2016 The MathWorks, Inc.

% Register variables in the workspace before the project is loaded
initVars = who;

% Variants Conditions
Variants.Command = 0;       % 0: Signal builder, 1: Joystick, 2: Pre-saved data, 3: Pre-saved data in a Spreadsheet
Variants.Sensors = 0;       % 0: Feedthrough, 1: Dynamics
Variants.Environment = 1;   % 0: Constant, 1: Variable
Variants.Vehicle = 1;       % 0: Linear dynamics, 1: Nonlinear dynamics
Variants.Visualization = 3; % 0: Scopes, 1: Send values to workspace, 2: FlightGear, 3: Simulink 3D.
Variants.Actuators = 0;     % 0: Feedthrough, 1: Linear second order, 2: Nonlinear 2nd order, 3: Linear first order.
Variants.Controller = 0;    % 0: Altitude hold, 1: Stabilize attitude for joystick.

% Add enum structure for the Variants
 Simulink.defineIntEnumType('Variants',{'Command','Vehicle','Environment',...
     'Sensors','Visualization','Actuators','Controller'},[0;0;0;0;0;0;0]);
 
% Bus definitions 
asbBusDefinitionCommand; 
asbBusDefinitionSensors;
asbBusDefinitionEnvironment;
asbBusDefinitionStates;

% Sampling rate
Ts= 0.01;

% Mass properties
mass = 0.078;
inertia = 1.0e-03*[0.3523 0.0001 0;... 
                   0.0001 0.3522 0;...
                   0      0 0.6788];

% Geometric properties
thrustArm = 0.10795;

% Initial conditions
initDate = [2015 1 1 0 0 0];
initPosLLA = [42.299886 -71.350447 71.3232];
initPosNED = [57 95 0];
initVb = [0 0 0];
initEuler = [0 0 0];
initAngRates = [0 0 0];

% Altitude hold value
altitudeHoldValue = initPosLLA(3) + 16.15;

%% Custom Variables
% Add your variables here:
% myvariable = 0;

% Register variables after the project is loaded and store the variables in
% initVars so they can be cleared later on the project shutdown.
endVars = who;
initVars = setdiff(endVars,initVars);
clear endVars;
