%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script loads mass and Inertia values defined by the user. Derivatives
%of the mass and Inertia tensor can be specified.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mass input
LD.Inertia.mass = 1;   %kg

% Gravity center vector position
LD.Inertia.CG = [LD.xcg_CMA , 0 , 0]'; 

% Inertia tensor matrix
LD.Inertia.Ixx = 1;      %kg*m^2
LD.Inertia.Iyy = 1;      %kg*m^2
LD.Inertia.Izz = 1;      %kg*m^2
LD.Inertia.Ixy = 0;      %kg*m^2
LD.Inertia.Ixz = 0;      %kg*m^2
LD.Inertia.Iyz = 0;      %kg*m^2

% Mass flow: Contains one or more rates of change of mass (positive if accreted, negative if ablated).
LD.Inertia.dmass = 0;    %kg/s

% CG position law function of mass. Defining the parameter LD.Inertia.mass,
% CG position is interpolated between these values:
LD.Inertia.fullMass    = LD.Inertia.mass;
LD.Inertia.emptyMass   = LD.Inertia.mass;
LD.Inertia.fullCG      = LD.Inertia.CG;
LD.Inertia.emptyCG     = LD.Inertia.CG;

% Contains the rate of change of Inertia tensor matrix.
LD.Inertia.dIxx = 0;     %kg*m^2/s
LD.Inertia.dIyy = 0;     %kg*m^2/s
LD.Inertia.dIzz = 0;     %kg*m^2/s
LD.Inertia.dIxy = 0;     %kg*m^2/s
LD.Inertia.dIxz = 0;     %kg*m^2/s
LD.Inertia.dIyz = 0;     %kg*m^2/s

% Contains one or more relative velocities at which the mass is accreted to
%or ablated from the body in body-fixed axes.
LD.Inertia.dmassdxb = 0; %m/s  Velocity component in x_body axis
LD.Inertia.dmassdyb = 0; %m/s  Velocity component in y_body axis
LD.Inertia.dmassdzb = 0; %m/s  Velocity component in z_body axis
