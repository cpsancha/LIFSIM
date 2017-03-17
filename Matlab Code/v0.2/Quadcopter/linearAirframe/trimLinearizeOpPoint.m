function trimLinearizeOpPoint()
% trimLinearizeOpPoint - Trim and Linearize aircraft model
% This function creates a linear system object that contains a linear model
% of the aircraft based on an trimmed operating point. The linear system
% object is assigned to linsys and saved to a .mat file. It uses Simulink
% Control Design(TM).

% Copyright 2013 The MathWorks, Inc.

%% Variants Conditions
oldVariants = evalin('base','Variants');

% Set up variants for trimming
Variants.Vehicle = 1;
Variants.Actuators = 3;

%Send variants to the workspace
assignin('base','Variants',Variants);

%% Obtain inputs and outputs from the model
model = 'trimNonlinearAirframe';
if ~bdIsLoaded(model)
    load_system(model);
end
io = getlinio(model);

%% Create the operating point specification object.
opspec = operspec(model);

%% Set the constraints on the states in the model.
% - The defaults for all states are Known = false, SteadyState = true,
%   Min = -Inf, and Max = Inf.

% State (1) - trimNonlinearAirframe/Airframe|nonlinearAirframe/Nonlinear/6DOF (Quaternion)/Calculate DCM & Euler Angles/q0 q1 q2 q3
% - Default model initial conditions are used to initialize optimization.

% State (2) - trimNonlinearAirframe/Airframe|nonlinearAirframe/Nonlinear/6DOF (Quaternion)/p,q,r
% - Default model initial conditions are used to initialize optimization.
opspec.States(2).Known = [true;true;true];

% State (3) - trimNonlinearAirframe/Airframe|nonlinearAirframe/Nonlinear/6DOF (Quaternion)/ub,vb,wb
% - Default model initial conditions are used to initialize optimization.
opspec.States(3).Known = [true;true;true];

% State (4) - trimNonlinearAirframe/Airframe|nonlinearAirframe/Nonlinear/6DOF (Quaternion)/xe,ye,ze
% - Default model initial conditions are used to initialize optimization.

% State (5) - trimNonlinearAirframe/Airframe|nonlinearAirframe/Nonlinear/AC model/Motor Forces and Torques/Thrust Dynamics/Linear 1st Order/Thruster 1
% - Default model initial conditions are used to initialize optimization.

% State (6) - trimNonlinearAirframe/Airframe|nonlinearAirframe/Nonlinear/AC model/Motor Forces and Torques/Thrust Dynamics/Linear 1st Order/Thruster 2
% - Default model initial conditions are used to initialize optimization.

% State (7) - trimNonlinearAirframe/Airframe|nonlinearAirframe/Nonlinear/AC model/Motor Forces and Torques/Thrust Dynamics/Linear 1st Order/Thruster 3
% - Default model initial conditions are used to initialize optimization.

% State (8) - trimNonlinearAirframe/Airframe|nonlinearAirframe/Nonlinear/AC model/Motor Forces and Torques/Thrust Dynamics/Linear 1st Order/Thruster 4
% - Default model initial conditions are used to initialize optimization.

% State (9) - trimNonlinearAirframe/Airframe|nonlinearAirframe/Nonlinear/AC model/Motor Forces and Torques/Torque Dynamics/Linear 1st Order/Thruster 1
% - Default model initial conditions are used to initialize optimization.

% State (10) - trimNonlinearAirframe/Airframe|nonlinearAirframe/Nonlinear/AC model/Motor Forces and Torques/Torque Dynamics/Linear 1st Order/Thruster 2
% - Default model initial conditions are used to initialize optimization.

% State (11) - trimNonlinearAirframe/Airframe|nonlinearAirframe/Nonlinear/AC model/Motor Forces and Torques/Torque Dynamics/Linear 1st Order/Thruster 3
% - Default model initial conditions are used to initialize optimization.

% State (12) - trimNonlinearAirframe/Airframe|nonlinearAirframe/Nonlinear/AC model/Motor Forces and Torques/Torque Dynamics/Linear 1st Order/Thruster 4
% - Default model initial conditions are used to initialize optimization.

%% Set the constraints on the inputs in the model.
% - The defaults for all inputs are Known = false, Min = -Inf, and
% Max = Inf.

% Input (1) - trimNonlinearAirframe/Actuators
% - Default model initial conditions are used to initialize optimization.
opspec.Inputs(1).Min = [1160;1160;1160;1160];
opspec.Inputs(1).Max = [1796;1796;1796;1796];

% Input (2) - trimNonlinearAirframe/Gravity ned
opspec.Inputs(2).u = [0;0;9.81];
opspec.Inputs(2).Known = [true;true;true];

% Input (3) - trimNonlinearAirframe/Air Temp
opspec.Inputs(3).u = 288;
opspec.Inputs(3).Known = true;

% Input (4) - trimNonlinearAirframe/speed sound
opspec.Inputs(4).u = 340;
opspec.Inputs(4).Known = true;

% Input (5) - trimNonlinearAirframe/pressure
opspec.Inputs(5).u = 101300;
opspec.Inputs(5).Known = true;

% Input (6) - trimNonlinearAirframe/air density
opspec.Inputs(6).u = 1.22;
opspec.Inputs(6).Known = true;

% Input (7) - trimNonlinearAirframe/Magnetic Field
% - Default model initial conditions are used to initialize optimization.
opspec.Inputs(7).Known = [true;true;true];

%% Set the constraints on the outputs in the model.
% - The defaults for all outputs are Known = false, Min = -Inf, and
% Max = Inf.

% Output (1) - trimNonlinearAirframe/V_body
% - Default model initial conditions are used to initialize optimization.
opspec.Outputs(1).Known = [true;true;true];

% Output (2) - trimNonlinearAirframe/Omega_body
% - Default model initial conditions are used to initialize optimization.
opspec.Outputs(2).Known = [true;true;true];

% Output (3) - trimNonlinearAirframe/Bus Selector
% - Default model initial conditions are used to initialize optimization.

% Output (4) - trimNonlinearAirframe/Bus Selector
% - Default model initial conditions are used to initialize optimization.

%% Create the options
opt = findopOptions('DisplayReport','iter');

%% Perform the operating point search.
[op,opreport] = findop(model,opspec,opt); %#okNASGU

% Linearize the model for the given operating and inputs/outputs.
linsys = linearize(model,io,op); %#okNASGU

% Restore Variants
assignin('base','Variants',oldVariants);

% Save trim points and linear model
p = simulinkproject;
save(fullfile(p.RootFolder,'linearAirframe','linearizedAirframe.mat'),...
    'linsys','op','opreport');

close_system(model);