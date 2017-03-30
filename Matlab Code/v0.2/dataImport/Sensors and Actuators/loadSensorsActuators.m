

%% Establish if sensors have bias and noise or not
% noisySensors = 0 --> The sensors do not have noise
% noisySensors = 1 --> The sensors have noise

    noisySensors = 1;
    
    
    
%% Establish the Actuators Dynamics
% actuatorsDynamics = 0 --> Feedthrough
% actuatorsDynamics = 1 --> Linear 1st order
% actuatorsDynamics = 2 --> Linear 2nd order
% actuatorsDynamics = 3 --> Nonlinear 2nd order

    actuatorsDynamics = 0;