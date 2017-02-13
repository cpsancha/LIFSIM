%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script loads mass and inertia values defined by the user. Derivatives
%of the mass and inertia tensor can be specified.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mass input
LD.inertia.mass = 1;   %kg

% Inertia tensor matrix
LD.inertia.Ixx = 1;      %kg*m^2
LD.inertia.Iyy = 1;      %kg*m^2
LD.inertia.Iyy = 1;      %kg*m^2
LD.inertia.Ixy = 0;      %kg*m^2
LD.inertia.Ixz = 0;      %kg*m^2
LD.inertia.Iyz = 0;      %kg*m^2

% Mass flow: Contains one or more rates of change of mass (positive if accreted, negative if ablated).
LD.inertia.dmass = 0;    %kg/s

% Contains the rate of change of inertia tensor matrix.
LD.inertia.dIxx = 0;     %kg*m^2/s
LD.inertia.dIyy = 0;     %kg*m^2/s
LD.inertia.dIyy = 0;     %kg*m^2/s
LD.inertia.dIxy = 0;     %kg*m^2/s
LD.inertia.dIxz = 0;     %kg*m^2/s
LD.inertia.dIyz = 0;     %kg*m^2/s

% Contains one or more relative velocities at which the mass is accreted to
%or ablated from the body in body-fixed axes.
LD.inertia.dmassdxb = 0; %m/s  Velocity component in x_body axis
LD.inertia.dmassdyb = 0; %m/s  Velocity component in y_body axis
LD.inertia.dmassdzb = 0; %m/s  Velocity component in z_body axis
