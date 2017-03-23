function [Voltage,Height,Vflight,Thrust,Power,Current,RPM,exitFlag] = propulsionDataMapping(motorData,...
          propellerData,voltageLim,heightLim,vflightLim,voltageLength,...
          heightLength,vflightLength,x0)
      
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
    exitFlag= zeros(voltageLength,heightLength,vflightLength);

%Define initial conditions
    if ~exist('x0','var')
        % x0 parameter does not exist, so default it to something
        x0 = [motorData.I0,7500,0,0];
    end


%Define @fsolve options
    options = optimoptions('fsolve',...
                           'Display','final-detailed',...
                           'TolFun',1e-6,...
                           'MaxIterations',1e3);
    
%Fill output variables
    parfor i=1:voltageLength
        for j=1:heightLength
            for k=1:vflightLength
                [X,~,exitflag,~] = fsolve(@(x)motorPropellerCoupling(x,Voltage(i),...
                    Height(j),Vflight(k),motorData,propellerData),x0,options); %#ok<PFBNS>
                exitFlag(i,j,k) = exitflag;
                Current(i,j,k) = X(1);
                RPM    (i,j,k) = X(2);
                Power  (i,j,k) = X(3);
                Thrust (i,j,k) = X(4);
                switch exitflag
                    case 2
                        disp(' ')
                        disp(['For values: Voltage=',num2str(Voltage(i)),' Height=',num2str(Height(j)),' Vflight=',num2str(Vflight(k))])
                        disp('Equation solved. Change in x smaller than the specified tolerance.')
                    case 3
                        disp(' ')
                        disp(['For values: Voltage=',num2str(Voltage(i)),' Height=',num2str(Height(j)),' Vflight=',num2str(Vflight(k))])
                        disp('Equation solved. Change in residual smaller than the specified tolerance.')
                    case 4
                        disp(' ')
                        disp(['For values: Voltage=',num2str(Voltage(i)),' Height=',num2str(Height(j)),' Vflight=',num2str(Vflight(k))])
                        disp('Equation solved. Magnitude of search direction smaller than specified tolerance.')
                    case 0
                        disp(' ')
                        disp(['For values: Voltage=',num2str(Voltage(i)),' Height=',num2str(Height(j)),' Vflight=',num2str(Vflight(k))])
                        disp('Number of iterations exceeded options.MaxIterations or number of function evaluations exceeded options.MaxFunctionEvaluations.')
                    case -1
                        disp(' ')
                        disp(['For values: Voltage=',num2str(Voltage(i)),' Height=',num2str(Height(j)),' Vflight=',num2str(Vflight(k))])
                        disp('Output function or plot function stopped the algorithm.')
                    case -2
                        disp(' ')
                        disp(['For values: Voltage=',num2str(Voltage(i)),' Height=',num2str(Height(j)),' Vflight=',num2str(Vflight(k))])
                        disp('Equation not solved. The exit message can have more information.')
                    case -3
                        disp(' ')
                        disp(['For values: Voltage=',num2str(Voltage(i)),' Height=',num2str(Height(j)),' Vflight=',num2str(Vflight(k))])
                        disp('Equation not solved. Trust region radius became too small (trust-region-dogleg algorithm).')
                end
            end
        end
    end

end

