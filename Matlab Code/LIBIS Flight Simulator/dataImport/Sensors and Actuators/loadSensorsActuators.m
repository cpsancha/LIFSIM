%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                                                         %
%   Script to load the pertinent data of the sensors and the actuators    %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Establish if sensors have bias and noise or not
% noisySensors = 0 --> The sensors do not have noise
% noisySensors = 1 --> The sensors have noise

    noisySensors = 1;
    
    
    
%% Establish the Actuators Dynamics
% actuatorsDynamics = 0 --> Feedthrough
% actuatorsDynamics = 1 --> Linear 1st order (LFO)
% actuatorsDynamics = 2 --> Linear 2nd order (LSO)
% actuatorsDynamics = 3 --> Nonlinear 2nd order (NLSO)

    actuatorsDynamics = 0;
    
  
    
    
%% Define the actuators transfer functions parameters
% Linear First Order Actuator
    LD.Actuators.LFO.timeLag = 0.05;
    
% Linear Second Order Actuator
    LD.Actuators.LSO.naturalFreq  = 44;
    LD.Actuators.LSO.dampingRatio = 0.707106781186547;
    LD.Actuators.LSO.initialPosition = 0;
    LD.Actuators.LSO.initialVelocity = 0;

% NonLinear Second Order Actuator
    LD.Actuators.NLSO.naturalFreq   = 0;
    LD.Actuators.NLSO.dampingRatio  = 0;
    LD.Actuators.NLSO.rateLimit     = 0;
    LD.Actuators.NLSO.maxDeflection = 0;
    LD.Actuators.NLSO.minDeflection = 0;
    LD.Actuators.NLSO.initialPosition = 0;
    LD.Actuators.NLSO.initialVelocity = 0;
    
    
    
    