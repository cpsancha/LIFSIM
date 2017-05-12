function [validData,Voltage,Height,Vflight,Thrust,Power,Current,RPM,exitFlag] = ...
    loadMappedData(filename,expectedMotorModel,expectedPropellerModel,...
    expectedVoltageLim,expectedHeightLim,expectedVflightLim,expectedVoltageLength,...
    expectedHeightLength,expectedVflightLength)
%LOADMAPPEDDATA Summary of this function goes here
    
    %load data contained in file
    load(filename);
    
    if isequal(motorModel,     expectedMotorModel)     && ...
       isequal(propellerModel, expectedPropellerModel) && ...
       isequal(voltageLim,     expectedVoltageLim)     && ...
       isequal(heightLim,      expectedHeightLim)      && ...
       isequal(vflightLim,     expectedVflightLim)     && ...
       isequal(voltageLength,  expectedVoltageLength)  && ...
       isequal(heightLength,   expectedHeightLength)   && ...
       isequal(vflightLength,  expectedVflightLength)
   
        validData = 1;
        Voltage = linspace(voltageLim(1),voltageLim(2),voltageLength);
        Height  = linspace( heightLim(1), heightLim(2), heightLength);
        Vflight = linspace(vflightLim(1),vflightLim(2),vflightLength);
        
    else
        validData = 0;
        Voltage = [];
        Height  = [];
        Vflight = [];
        Thrust  = [];
        Power   = [];
        Current = [];
        RPM     = [];
        exitFlag= [];
    end

end

