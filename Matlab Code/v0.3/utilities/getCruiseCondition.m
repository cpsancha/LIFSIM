VEAS=23.11;
height =91.44;
[~,~,~,rho] =atmosisa(height);



Cl=2*LD.Inertia.mass*9.8/(rho*LD.sref*VEAS^2);

%% Control variables initialization
deltae_degrees=0;

deltat_2=0;  % motores delanteros
deltat_3=0;  % motores traseros
deltat_4=0;  % motores delanteros
deltat_5=0;  % motores traseros
deltat_1=1;   % ala fija

%% Interpolación de las derivadas de estabilidad para la condición especificada
X1 = LD.alpha;
X2 = LD.beta;
X3 = LD.alt;
X4 = LD.xcg;
X5 = LD.deltae;
X6 = LD.
Xq1=-0*pi/180;
Xq2=0;
V=LD.Stability.CL0(:,:,1,1,1,1,1,1);


    %Array of altitudes.
        LD.alt    = 91.44;
    %Array of positions of the center of gravity. Distance from the nose of
    %the plane to center of gravity (meters) in body axis with origin at the
    %nose
        LD.xcg    = -1.013;
    %Array of elevator, rudder, and right and left flaperons streamwise control 
    %deflection angles, which are defined positive for trailing-edge down.
        LD.deltae  = 10*pi/180;
        LD.deltar  = 10*pi/180;
        LD.deltafr = 10*pi/180;
        LD.deltafl = 10*pi/180;

clcero=interpn(X1,X2,V,Xq1,Xq2)

Cl0=LD.Stability.CL0(1,1,1,1,1,1,1);
Cl_alpha=LD.Stability.CLalpha(1,1,1,1,1,1,1);
Cl_alpha_dot=LD.Stability.CLalpha_dot(1,1,1,1,1,1,1);
Cm_alpha=LD.Stability.Cmalpha(1,1,1,1,1,1,1);
Cl_delta_e=LD.Stability.CLdeltae(1,1,1,1,1,1,1);
Cm_delta_e=LD.Stability.Cmdeltae(1,1,1,1,1,1,1);
Cm0=LD.Stability.Cm0(1,1,1,1,1,1,1);

alpha_wb=(Cl-Cl0+Cl_delta_e*Cm0/Cm_delta_e)./(Cl_alpha-Cl_delta_e*Cm_alpha/Cm_delta_e);
delta_e=-(Cm0+Cm_alpha*alpha_wb)/Cm_delta_e;


% alphacruise=(CL-CL0)/CLalpha*180/pi;

% 

% initialValues.u0=cos(alpha_wb)*VEAS;
% initialValues.w0=sin(alpha_wb)*VEAS;


initialValues.Ze0 = -91.44; %Cruise height
initialValues.u0 = 23.1;     %Cruise speeed

initialValues.pitch0 = alpha_wb;

%Forces
% alpha_wb=0.171;
qbar=0.5*rho*VEAS^2;
CL=(Cl0+alpha_wb*Cl_alpha+Cl_delta_e*delta_e);
L=qbar*LD.sref*CL;
