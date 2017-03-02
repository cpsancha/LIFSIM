function [ Error ] = motorPropellerCoupling(x,Vesc,h,Vflight,motorData,propellerData)
%MOTORPROPELLERCOUPLING Summary of this function goes here:
%   INPUTS:
%    x = Array of unknown values
%     x(1) = Current of the motor in Amps (Im)
%     x(2) = RPM of the motor-propeller in rev/min (RPM)
%     x(3) = Power delivered by the motor to the shaft in Watts (Pmotor)
%     x(4) = Power absorbed by the propeller in Watts (Pprop)
%     x(5) = Thrust provided by the propeller in Newtons (T)
%    Vesc = Voltage from the ESC (Electronic Speed Control) in Volts
%    h    = Height of flight in meters
%    Vfl  = Flight speed in m/s
%    motorData = Structure of data from the motor
%    propellerData = Structure of data from the propeller  
%   OUTPUTS:
%    Error = Solution to make zero


%% Variable rename for better understanding
    Im  = x(1); 
    RPM = x(2);
    P   = x(3);
    T   = x(4);
    RPS = RPM/60;
    [~,~,~,rho] = atmosisa(h);
    [VFLIGHT,RPMmesh] = meshgrid(propellerData.V,propellerData.RPM);
    Ct = interp2(VFLIGHT,RPMmesh,propellerData.Ct,Vflight,RPM);
    Cp = interp2(VFLIGHT,RPMmesh,propellerData.Cp,Vflight,RPM);
    
%% Motor Eq:
    Error(1) = Im - motorData.I0 - (P/(RPM/motorData.Kv));
    Error(2) = Vesc - (RPM/motorData.Kv) - (Im*motorData.Rm);

%% Propeller Eq:
    Error(3) = P - rho*RPS^3*propellerData.Diameter^5*Cp;
    Error(4) = T - rho*RPS^2*propellerData.Diameter^4*Ct;

    
end

