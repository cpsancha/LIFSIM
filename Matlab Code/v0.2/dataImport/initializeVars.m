%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to initialize all the necessary variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run('loadData.m')

% run('loadBusData')

%Direction of flat Earth x-axis (degrees clockwise from north):
    initialValues.heading0 = 0;
    
%Initial geodetic latitude and longitude [deg]:
    initialValues.LatLong0 = [45 120];
    
%Reference height from surface of Earth to flat Earth frame with regard to 
%Earth frame, in same units as flat Earth position (meters):
    initialValues.href = 0;
    
%Initial position in inertial axes [Xe,Ye,Ze]:
    initialValues.Xe0 = 0;
    initialValues.Ye0 = 0;
    initialValues.Ze0 = -0.00391527660175;
    
%Initial velocity in body axes [u,v,w]:
    initialValues.u0 = 0;
    initialValues.v0 = 0;
    initialValues.w0 = 0;
%Initial Euler orientation [roll, pitch, yaw]:
    initialValues.roll0  = 0;
    initialValues.pitch0 = 0;
    initialValues.yaw0   = 0;
%Initial body rotation rates [p,q,r]:
    initialValues.p0 = 0;
    initialValues.q0 = 0;
    initialValues.r0 = 0;
    
%% Ground reaction constants
    groundReaction.K = 1e6;
    groundReaction.D = 2000;
    groundReaction.Mu = 0;%0.1;
    
%% Control variables initialization

    deltae_degrees  = 0;
    deltar_degrees  = 0;
    deltafr_degrees = 0;
    deltafl_degrees = 0;
    delta_e=0;

    deltat_2=0;  % motores delanteros
    deltat_3=0;  % motores traseros
    deltat_4=0;  % motores delanteros
    deltat_5=0;  % motores traseros
    deltat_1=1;   % ala fija