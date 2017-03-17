function asbBusDefinitionSensors() 
% ASBBUSDEFINITIONSENSORS initializes a set of bus objects in the MATLAB base workspace 
% Copyright 2013-2015 The MathWorks, Inc.

% Bus object: Sensors 
clear elems;
elems(1) = Simulink.BusElement;
elems(1).Name = 'AccelMeas_body';
elems(1).Dimensions = [3 1];
elems(1).DimensionsMode = 'Fixed';
elems(1).DataType = 'double';
elems(1).SampleTime = -1;
elems(1).Complexity = 'real';
elems(1).SamplingMode = 'Sample based';
elems(1).Min = [];
elems(1).Max = [];
elems(1).DocUnits = '';
elems(1).Description = '';

elems(2) = Simulink.BusElement;
elems(2).Name = 'OmegaMeas_body';
elems(2).Dimensions = [3 1];
elems(2).DimensionsMode = 'Fixed';
elems(2).DataType = 'double';
elems(2).SampleTime = -1;
elems(2).Complexity = 'real';
elems(2).SamplingMode = 'Sample based';
elems(2).Min = [];
elems(2).Max = [];
elems(2).DocUnits = '';
elems(2).Description = '';

elems(3) = Simulink.BusElement;
elems(3).Name = 'LLA';
elems(3).Dimensions = [1 3];
elems(3).DimensionsMode = 'Fixed';
elems(3).DataType = 'double';
elems(3).SampleTime = -1;
elems(3).Complexity = 'real';
elems(3).SamplingMode = 'Sample based';
elems(3).Min = [];
elems(3).Max = [];
elems(3).DocUnits = '';
elems(3).Description = '';

SensorsBus = Simulink.Bus;
SensorsBus.HeaderFile = '';
SensorsBus.Description = '';
SensorsBus.DataScope = 'Auto';
SensorsBus.Alignment = -1;
SensorsBus.Elements = elems;
assignin('base','SensorsBus', SensorsBus)

