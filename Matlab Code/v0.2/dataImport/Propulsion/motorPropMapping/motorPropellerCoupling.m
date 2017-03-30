function [ Error ] = motorPropellerCoupling(x,Vesc,h,Vflight,motorData,propellerData)
%MOTORPROPELLERCOUPLING Summary of this function goes here:
%   INPUTS:
%    x = Array of unknown values
%     x(1).*1e1 = Current of the motor in Amps (Im)
%     x(2).*1e4 = RPM of the motor-propeller in rev/min (RPM)
%     x(3).*1e1 = Power delivered by the motor to the shaft in Watts (Pmotor)
%     x(3).*1e1 = Power absorbed by the propeller in Watts (Pprop)
%     x(4).*1e0 = Thrust provided by the propeller in Newtons (T)
%    Vesc = Voltage from the ESC (Electronic Speed Control) in Volts
%    h    = Height of flight in meters
%    Vfl  = Flight speed in m/s
%    motorData = Structure of data from the motor
%    propellerData = Structure of data from the propeller  
%   OUTPUTS:
%    Error = Solution to make zero


%% Variable rename for better understanding
    Im  = x(1).*1e1; 
    RPM = x(2).*1e4;
    P   = x(3).*1e2;
    T   = x(4).*1e1;
    RPS = RPM./60;
    [~,~,~,rho] = atmosisa(h);
    [VFLIGHT,RPMmesh] = meshgrid(propellerData.V,propellerData.RPM);
    Ct = interp2(VFLIGHT,RPMmesh,propellerData.Ct,Vflight,RPM,'linear',0);
    Cp = interp2(VFLIGHT,RPMmesh,propellerData.Cp,Vflight,RPM,'linear',0);

    
%% Solve Nonlinear System
    if Vesc >= motorData.I0*motorData.Rm
        %Motor Eq:
            Error(1) = (Im - motorData.I0) - ((P.*motorData.Kv)./RPM);
            Error(2) = Vesc - (RPM./motorData.Kv) - (Im.*motorData.Rm);
        %Propeller Eq:
            Error(3) = P - rho.*RPS.^3.*propellerData.Diameter.^5.*Cp;
            Error(4) = T - rho.*RPS.^2.*propellerData.Diameter.^4.*Ct;
    else
        %Motor Eq:
            Error(1) = Im;
            Error(2) = RPM;
        %Propeller Eq:
            Error(3) = P;
            Error(4) = T;  
    end
    
end

