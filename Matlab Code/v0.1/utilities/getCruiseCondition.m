VEAS=23.11;
rho=1.28;
Cl=2*LD.Inertia.mass*9.8/(rho*LD.sref*VEAS^2);

Cl0=LD.Stability.CL0(1,1,1,1,1,1,1);
Cl_alpha=LD.Stability.CLalpha(1,1,1,1,1,1,1);
Cl_alpha_dot=LD.Stability.CLalpha_dot(1,1,1,1,1,1,1);
Cm_alpha=LD.Stability.Cmalpha(1,1,1,1,1,1,1);
Cl_delta_e=LD.Stability.CLdeltae(1,1,1,1,1,1,1);
Cm_delta_e=LD.Stability.Cmdeltae(1,1,1,1,1,1,1);
Cm0=LD.Stability.Cm0(1,1,1,1,1,1,1);

alpha_wb=(Cl-Cl0+Cl_delta_e*Cm0/Cm_delta_e)./(Cl_alpha-Cl_delta_e*Cm_alpha/Cm_delta_e);
delta_e=-(Cm0+Cm_alpha*alpha_wb)/Cm_delta_e;
delta_e=1.2*delta_e;

% alphacruise=(CL-CL0)/CLalpha*180/pi;

% 

initialValues.u0=cos(alpha_wb)*VEAS;
initialValues.w0=sin(alpha_wb)*VEAS;


initialValues.Ze0 = -91.44; %Cruise height

initialValues.pitch0 = alpha_wb;

%Forces
% alpha_wb=0.171;
qbar=0.5*rho*VEAS^2;
CL=(Cl0+alpha_wb*Cl_alpha+Cl_delta_e*delta_e);
L=qbar*LD.sref*CL;
