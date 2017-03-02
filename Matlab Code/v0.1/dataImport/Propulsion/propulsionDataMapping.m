function [Thrust,Power,Current,RPM] = propulsionDataMapping(motorData,...
          propellerData,voltageLim,heightLim,vflightLim,voltageLength,...
          heightLength,vflightLength)
      
%PROPULSIONDATAMAPPING Summary of this function goes here
%% MAP THRUST, POWER, CURRENT AND RPM IN FUNCTION OF VOLTAGE, HEIGHT AND VELOCITY
%Vectors:    
    Voltage = linspace(voltageLim(1),voltageLim(2),voltageLength);
    Height  = linspace( heightLim(1), heightLim(2), heightLength);
    Vflight = linspace(vflightLim(1),vflightLim(2),vflightLength);

%Define output variables:
    Thrust  = zeros(voltageLength,heightLength,vflightLength);
    Power   = zeros(voltageLength,heightLength,vflightLength);
    Current = zeros(voltageLength,heightLength,vflightLength);
    RPM     = zeros(voltageLength,heightLength,vflightLength);

%Define initial conditions
    x0 = [20,10000,200,100];
    
%Fill output variables
    for i=1:voltageLength
       for j=1:heightLength
          for k=1:vflightLength
              X = fsolve(@(x)motorPropellerCoupling(x,Voltage(i),Height(j),...
                         Vflight(k),motorData,propellerData),x0);
              Current(i,j,k) = X(1);
              RPM    (i,j,k) = X(2);
              Power  (i,j,k) = X(3);
              Thrust (i,j,k) = X(4);
          end
       end
    end

end

