%% Vehicle Nonlinear Variables
% Copyright 2013-2015 The MathWorks, Inc.

% 6DOF
Vehicle.Nonlinear.Nonlinear.SixDOF.initGreenwich = 0;
Vehicle.Nonlinear.Nonlinear.SixDOF.quatGain = 1;
% Actuators
% Linear
Vehicle.Nonlinear.Nonlinear.ACModel.Actuator.Linear.natFreq = 1;
Vehicle.Nonlinear.Nonlinear.ACModel.Actuator.Linear.damping = 0.3;
Vehicle.Nonlinear.Nonlinear.ACModel.Actuator.Linear.initPosition = 0;
Vehicle.Nonlinear.Nonlinear.ACModel.Actuator.Linear.initVelocity = 0;
% Nonlinear
Vehicle.Nonlinear.Nonlinear.ACModel.Actuator.NonLinear.natFreq = 1;
Vehicle.Nonlinear.Nonlinear.ACModel.Actuator.NonLinear.damping = 0.3;
Vehicle.Nonlinear.Nonlinear.ACModel.Actuator.NonLinear.maxDeflection = 20*pi/180;
Vehicle.Nonlinear.Nonlinear.ACModel.Actuator.NonLinear.minDeflection = 20*pi/180;
Vehicle.Nonlinear.Nonlinear.ACModel.Actuator.NonLinear.rateLimit = 500*pi/180;
Vehicle.Nonlinear.Nonlinear.ACModel.Actuator.NonLinear.initPosition = 1;
Vehicle.Nonlinear.Nonlinear.ACModel.Actuator.NonLinear.initVelocity = 0;
% Position on Earth
Vehicle.Nonlinear.Nonlinear.PositionOnEarth.href = -71.3232;
Vehicle.Nonlinear.Nonlinear.PositionOnEarth.FlatEarthToLLA.xAxis = 90;
% First Order Linear
Vehicle.Nonlinear.Nonlinear.ACModel.Actuator.LinearFirst.invTs = 15;
% Look-up Tables
% Thrust
Vehicle.Nonlinear.Nonlinear.ACModel.ForcesMoments.Thrust.tableData = 0.00980665*[0 115 220 306 419 512 613 719 822];
Vehicle.Nonlinear.Nonlinear.ACModel.ForcesMoments.Thrust.breakpoints = [1160 1240 1339 1429 1528 1594 1660 1720 1796];
% Torque
Vehicle.Nonlinear.Nonlinear.ACModel.ForcesMoments.Torque.tableData = 0.00980665*0.3048*[0 10 20 35 42 58 63];
Vehicle.Nonlinear.Nonlinear.ACModel.ForcesMoments.Torque.breakpoints = [1160 1300 1389 1500 1600 1720 1796];
% Saturation
Vehicle.Nonlinear.Nonlinear.ACModel.ForcesMoments.Torque.lowerSat = 0;
Vehicle.Nonlinear.Nonlinear.ACModel.ForcesMoments.Torque.upperSat = 0.189;
Vehicle.Nonlinear.Nonlinear.ACModel.ForcesMoments.Thrust.lowerSat = 0;
Vehicle.Nonlinear.Nonlinear.ACModel.ForcesMoments.Thrust.upperSat = 8.1;
% Drag Calculation
Vehicle.Nonlinear.Nonlinear.ACModel.dragCalc.dragCoeff = [.1 .1 .3];
Vehicle.Nonlinear.Nonlinear.ACModel.dragCalc.diameter = 0.060293;