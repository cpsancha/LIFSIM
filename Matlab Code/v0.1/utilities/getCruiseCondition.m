VEAS=23.11;
Cl=2*LD.Inertia.mass*9.8/(1.225*LD.sref*VEAS^2);

Cl0=LD.Stability.CL0(1,1,1,1,1,1,1);
Cl_alpha=LD.Stability.CLalpha(1,1,1,1,1,1,1);
Cm_alpha=LD.Stability.Cmalpha(1,1,1,1,1,1,1);
Cl_delta_e=LD.Stability.CLdeltae(1,1,1,1,1,1,1);
Cm_delta_e=LD.Stability.Cmdeltae(1,1,1,1,1,1,1);
Cm0=LD.Stability.Cm0(1,1,1,1,1,1,1);

alpha_wb=(Cl-Cl0+Cl_delta_e*Cm0/Cm_delta_e)./(Cl_alpha-Cl_delta_e*Cm_alpha/Cm_delta_e);
delta_e=-(Cm0+Cm_alpha*alpha_wb)/Cm_delta_e;

% alphacruise=(CL-CL0)/CLalpha*180/pi;


initialValues.u0=cos(alpha_wb)*VEAS;
initialValues.v0=-sin(alpha_wb)*VEAS;

initialValues.Ze0 = -91.44; %Cruise height
