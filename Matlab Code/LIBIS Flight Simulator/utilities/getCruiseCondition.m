%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  getCruiseConditions  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gets the trim variables neede to perform an horizontal, straight 
% and uniform fligth, according to the velocity and height selected in the
% section 'input'. This variables are then overwritten over the initial
% values of the simulation

%% Input:

VEAS=23.11;
height =30;

%%
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
X6 = LD.deltar;
X7 = LD.deltafr;
X8 = LD.deltafl;

% Xqn = Punto en el que queremos conocer el valor de las derivadas.
Xq1 = 0; %Se asume que el equilibirio se encuentra para ángulos de ataque muy pequeños (zona lineal)
Xq2 = 0;
Xq3 = height+91.44;  %Debido a que el interpn no permite el clip, los valores no pueden ser inferiores a la altura para la que se definen los ceoficientes
Xq4 = LD.Inertia.CG(1);
Xq5 = 10*pi/180;
Xq6 = 0;
Xq7 = 10*pi/180;
Xq8 = 10*pi/180;

V = LD.Stability.CL0;
Cl0 = interpn(X1,X2,X3,X4,X5,X6,X7,X8,V,Xq1,Xq2,Xq3,Xq4,Xq5,Xq6,Xq7,Xq8);

V = LD.Stability.CLalpha;
Cl_alpha = interpn(X1,X2,X3,X4,X5,X6,X7,X8,V,Xq1,Xq2,Xq3,Xq4,Xq5,Xq6,Xq7,Xq8);

V = LD.Stability.Cmalpha;
Cm_alpha = interpn(X1,X2,X3,X4,X5,X6,X7,X8,V,Xq1,Xq2,Xq3,Xq4,Xq5,Xq6,Xq7,Xq8);

V = LD.Stability.CLdeltae;
Cl_delta_e = interpn(X1,X2,X3,X4,X5,X6,X7,X8,V,Xq1,Xq2,Xq3,Xq4,Xq5,Xq6,Xq7,Xq8);

V = LD.Stability.Cmdeltae;
Cm_delta_e = interpn(X1,X2,X3,X4,X5,X6,X7,X8,V,Xq1,Xq2,Xq3,Xq4,Xq5,Xq6,Xq7,Xq8);

V = LD.Stability.Cm0;
Cm0 = interpn(X1,X2,X3,X4,X5,X6,X7,X8,V,Xq1,Xq2,Xq3,Xq4,Xq5,Xq6,Xq7,Xq8);

%Drag derivatives:
V = LD.Stability.CD0;
CD0 = interpn(X1,X2,X3,X4,X5,X6,X7,X8,V,Xq1,Xq2,Xq3,Xq4,Xq5,Xq6,Xq7,Xq8);

V = LD.Stability.CDalpha;
CDalpha = interpn(X1,X2,X3,X4,X5,X6,X7,X8,V,Xq1,Xq2,Xq3,Xq4,Xq5,Xq6,Xq7,Xq8);

V = LD.Stability.CDdeltae;
CDdeltae = interpn(X1,X2,X3,X4,X5,X6,X7,X8,V,Xq1,Xq2,Xq3,Xq4,Xq5,Xq6,Xq7,Xq8);


%% Se resuelve la ecuación de trimado de forma analítica

alpha_wb=(Cl-Cl0+Cl_delta_e*Cm0/Cm_delta_e)./(Cl_alpha-Cl_delta_e*Cm_alpha/Cm_delta_e);
delta_e=-(Cm0+Cm_alpha*alpha_wb)/Cm_delta_e;


%% Equilibrio de fuerzas

qbar = 0.5*rho*VEAS^2;
CL = (Cl0+alpha_wb*Cl_alpha+Cl_delta_e*delta_e);
L  = qbar*LD.sref*CL;

CD = 0.0241+0.065*CL; %Polar roskam
CD = CD0+alpha_wb*CDalpha+CDdeltae*delta_e; %Derivadas de estabilidad
D = qbar*LD.sref*CD; 

%% Interpolacion en el mapeado del motor

% Se obtiene el empuje en funcion del voltaje de entrada al motor para la
% condición de vuelo (VEAS,height)

X1 = LD.Propulsion.motor1.Voltage;
X2 = LD.Propulsion.motor1.Height;
X3 = LD.Propulsion.motor1.Vflight;

V = LD.Propulsion.motor1.Thrust;

Xq1 = LD.Propulsion.motor1.Voltage;
Xq2 = height;
Xq3 = VEAS*cos(alpha_wb);

T_v = interpn (X1,X2,X3,V,Xq1,Xq2,Xq3);

% Se vuelve a interpolar para encontrar el voltaje que corresponde a la
% tracción necesaria en la condición de vuelo (T=D)

 [T_v, index] = unique(T_v); 
 Voltage = interp1(T_v, LD.Propulsion.motor1.Voltage(index), D);

% X1  = T_v;
% V   = LD.Propulsion.motor1.Voltage;
% Xq1 = D;
% 
% Voltage = interpn(X1,V,Xq1);
% 
% [~ , index] = min(abs(T_v-D));
% Voltage =  LD.Propulsion.motor1.Voltage(index)

deltat_1 = Voltage/(LD.Propulsion.ESC1.eta*29.6*0.98);

%% Se sobreescriben los valores iniciales a la condición de crucero:

initialValues.u0 = cos(alpha_wb)*VEAS;
initialValues.w0 = sin(alpha_wb)*VEAS;

initialValues.uned0 = VEAS;

initialValues.Ze0 = -height; %Cruise height
initialValues.pitch0 = alpha_wb; %Para que el vuelo sea horizontal, además de uniforme

initialValues.alpha_wb = alpha_wb;
initialValues.D = D;
initialValues.L = L;

% initialValues.u0 = 0;%cos(alpha_wb)*VEAS;
% initialValues.w0 = 0; %sin(alpha_wb)*VEAS;
% 
% initialValues.uned0 = 0;VEAS;
% 
% initialValues.Ze0 = -0.5;%-height; %Cruise height
% initialValues.pitch0 = 0; %alpha_wb; %Para que el vuelo sea horizontal, además de uniforme