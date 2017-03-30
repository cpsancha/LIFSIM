
%% Set values for gains in FCS
LD.FCS.Quadcopter.Pitch.dPitch       = 1;
LD.FCS.Quadcopter.Roll.dRoll         = 1;
LD.FCS.Quadcopter.Yaw.dYaw           = 1;% -0.7;
LD.FCS.Quadcopter.Throttle.dThrottle = 1;%-0.1;
LD.FCS.Quadcopter.Throttle.bias      = 0;

% LD.FCS.Quadcopter.Throttle.dThrottle = -0.004;%-0.1;
% LD.FCS.Quadcopter.Throttle.bias  = 1;
% LD.FCS.Quadcopter.Pitch.dPitch =       0.001;
% LD.FCS.Quadcopter.Roll.dRoll   =       0.001;
% LD.FCS.Quadcopter.Yaw.dYaw     =      -6;% -0.7;

    
%% Choose the flight condition
% LD.FCS.flightCondition = 0 --> LIBIS is acting as quadcopter
% LD.FCS.flightCondition = 1 --> Transition from quadcopter to fixed-wing plane
% LD.FCS.flightCondition = 2 --> LIBIS is acting as fixed-wing plane
% LD.FCS.flightCondition = 3 --> Transition from fixed-wing plane to quadcopter
% LD.FCS.flightCondition = 4 --> Some random conditions to test TO BE DONE....
  LD.FCS.fligthCondition = 0;
