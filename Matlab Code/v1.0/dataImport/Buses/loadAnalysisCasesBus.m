function loadAnalysisCasesBus() 
% BUSESDEFINITION initializes a set of bus objects in the MATLAB base workspace 

% Bus object: analysisCasesBus 
clear elems;
elems(1) = Simulink.BusElement;
elems(1).Name = 'alpha';
elems(1).Dimensions = 1;
elems(1).DimensionsMode = 'Fixed';
elems(1).DataType = 'Bus: alphaBus';
elems(1).SampleTime = -1;
elems(1).Complexity = 'real';
elems(1).Min = [];
elems(1).Max = [];
elems(1).DocUnits = '';
elems(1).Description = '';

elems(2) = Simulink.BusElement;
elems(2).Name = 'beta';
elems(2).Dimensions = 1;
elems(2).DimensionsMode = 'Fixed';
elems(2).DataType = 'Bus: betaBus';
elems(2).SampleTime = -1;
elems(2).Complexity = 'real';
elems(2).Min = [];
elems(2).Max = [];
elems(2).DocUnits = '';
elems(2).Description = '';

elems(3) = Simulink.BusElement;
elems(3).Name = 'alt';
elems(3).Dimensions = 1;
elems(3).DimensionsMode = 'Fixed';
elems(3).DataType = 'Bus: altBus';
elems(3).SampleTime = -1;
elems(3).Complexity = 'real';
elems(3).Min = [];
elems(3).Max = [];
elems(3).DocUnits = '';
elems(3).Description = '';

elems(4) = Simulink.BusElement;
elems(4).Name = 'xcg';
elems(4).Dimensions = 1;
elems(4).DimensionsMode = 'Fixed';
elems(4).DataType = 'Bus: xcgBus';
elems(4).SampleTime = -1;
elems(4).Complexity = 'real';
elems(4).Min = [];
elems(4).Max = [];
elems(4).DocUnits = '';
elems(4).Description = '';

elems(5) = Simulink.BusElement;
elems(5).Name = 'deltae';
elems(5).Dimensions = 1;
elems(5).DimensionsMode = 'Fixed';
elems(5).DataType = 'Bus: deltaeBus';
elems(5).SampleTime = -1;
elems(5).Complexity = 'real';
elems(5).Min = [];
elems(5).Max = [];
elems(5).DocUnits = '';
elems(5).Description = '';

elems(6) = Simulink.BusElement;
elems(6).Name = 'deltar';
elems(6).Dimensions = 1;
elems(6).DimensionsMode = 'Fixed';
elems(6).DataType = 'Bus: deltarBus';
elems(6).SampleTime = -1;
elems(6).Complexity = 'real';
elems(6).Min = [];
elems(6).Max = [];
elems(6).DocUnits = '';
elems(6).Description = '';

elems(7) = Simulink.BusElement;
elems(7).Name = 'deltafr';
elems(7).Dimensions = 1;
elems(7).DimensionsMode = 'Fixed';
elems(7).DataType = 'Bus: deltafrBus';
elems(7).SampleTime = -1;
elems(7).Complexity = 'real';
elems(7).Min = [];
elems(7).Max = [];
elems(7).DocUnits = '';
elems(7).Description = '';

elems(8) = Simulink.BusElement;
elems(8).Name = 'deltafl';
elems(8).Dimensions = 1;
elems(8).DimensionsMode = 'Fixed';
elems(8).DataType = 'Bus: deltaflBus';
elems(8).SampleTime = -1;
elems(8).Complexity = 'real';
elems(8).Min = [];
elems(8).Max = [];
elems(8).DocUnits = '';
elems(8).Description = '';

analysisCasesBus = Simulink.Bus;
analysisCasesBus.HeaderFile = '';
analysisCasesBus.Description = '';
analysisCasesBus.DataScope = 'Auto';
analysisCasesBus.Alignment = -1;
analysisCasesBus.Elements = elems;
clear elems;
assignin('base','analysisCasesBus', analysisCasesBus);

% Bus object: alphaBus 
clear elems;
elems(1) = Simulink.BusElement;
elems(1).Name = 'idxalpha';
elems(1).Dimensions = 1;
elems(1).DimensionsMode = 'Fixed';
elems(1).DataType = 'int32';
elems(1).SampleTime = -1;
elems(1).Complexity = 'real';
elems(1).Min = [];
elems(1).Max = [];
elems(1).DocUnits = '';
elems(1).Description = '';

elems(2) = Simulink.BusElement;
elems(2).Name = 'falpha';
elems(2).Dimensions = 1;
elems(2).DimensionsMode = 'Fixed';
elems(2).DataType = 'double';
elems(2).SampleTime = -1;
elems(2).Complexity = 'real';
elems(2).Min = [];
elems(2).Max = [];
elems(2).DocUnits = '';
elems(2).Description = '';

alphaBus = Simulink.Bus;
alphaBus.HeaderFile = '';
alphaBus.Description = '';
alphaBus.DataScope = 'Auto';
alphaBus.Alignment = -1;
alphaBus.Elements = elems;
clear elems;
assignin('base','alphaBus', alphaBus);

% Bus object: betaBus 
clear elems;
elems(1) = Simulink.BusElement;
elems(1).Name = 'idxbeta';
elems(1).Dimensions = 1;
elems(1).DimensionsMode = 'Fixed';
elems(1).DataType = 'int32';
elems(1).SampleTime = -1;
elems(1).Complexity = 'real';
elems(1).Min = [];
elems(1).Max = [];
elems(1).DocUnits = '';
elems(1).Description = '';

elems(2) = Simulink.BusElement;
elems(2).Name = 'fbeta';
elems(2).Dimensions = 1;
elems(2).DimensionsMode = 'Fixed';
elems(2).DataType = 'double';
elems(2).SampleTime = -1;
elems(2).Complexity = 'real';
elems(2).Min = [];
elems(2).Max = [];
elems(2).DocUnits = '';
elems(2).Description = '';

betaBus = Simulink.Bus;
betaBus.HeaderFile = '';
betaBus.Description = '';
betaBus.DataScope = 'Auto';
betaBus.Alignment = -1;
betaBus.Elements = elems;
clear elems;
assignin('base','betaBus', betaBus);

% Bus object: altBus 
clear elems;
elems(1) = Simulink.BusElement;
elems(1).Name = 'idxalt';
elems(1).Dimensions = 1;
elems(1).DimensionsMode = 'Fixed';
elems(1).DataType = 'int32';
elems(1).SampleTime = -1;
elems(1).Complexity = 'real';
elems(1).Min = [];
elems(1).Max = [];
elems(1).DocUnits = '';
elems(1).Description = '';

elems(2) = Simulink.BusElement;
elems(2).Name = 'falt';
elems(2).Dimensions = 1;
elems(2).DimensionsMode = 'Fixed';
elems(2).DataType = 'double';
elems(2).SampleTime = -1;
elems(2).Complexity = 'real';
elems(2).Min = [];
elems(2).Max = [];
elems(2).DocUnits = '';
elems(2).Description = '';

altBus = Simulink.Bus;
altBus.HeaderFile = '';
altBus.Description = '';
altBus.DataScope = 'Auto';
altBus.Alignment = -1;
altBus.Elements = elems;
clear elems;
assignin('base','altBus', altBus);

% Bus object: xcgBus 
clear elems;
elems(1) = Simulink.BusElement;
elems(1).Name = 'idxxcg';
elems(1).Dimensions = 1;
elems(1).DimensionsMode = 'Fixed';
elems(1).DataType = 'int32';
elems(1).SampleTime = -1;
elems(1).Complexity = 'real';
elems(1).Min = [];
elems(1).Max = [];
elems(1).DocUnits = '';
elems(1).Description = '';

elems(2) = Simulink.BusElement;
elems(2).Name = 'fxcg';
elems(2).Dimensions = 1;
elems(2).DimensionsMode = 'Fixed';
elems(2).DataType = 'double';
elems(2).SampleTime = -1;
elems(2).Complexity = 'real';
elems(2).Min = [];
elems(2).Max = [];
elems(2).DocUnits = '';
elems(2).Description = '';

xcgBus = Simulink.Bus;
xcgBus.HeaderFile = '';
xcgBus.Description = '';
xcgBus.DataScope = 'Auto';
xcgBus.Alignment = -1;
xcgBus.Elements = elems;
clear elems;
assignin('base','xcgBus', xcgBus);

% Bus object: deltaeBus 
clear elems;
elems(1) = Simulink.BusElement;
elems(1).Name = 'idxdeltae';
elems(1).Dimensions = 1;
elems(1).DimensionsMode = 'Fixed';
elems(1).DataType = 'int32';
elems(1).SampleTime = -1;
elems(1).Complexity = 'real';
elems(1).Min = [];
elems(1).Max = [];
elems(1).DocUnits = '';
elems(1).Description = '';

elems(2) = Simulink.BusElement;
elems(2).Name = 'fdeltae';
elems(2).Dimensions = 1;
elems(2).DimensionsMode = 'Fixed';
elems(2).DataType = 'double';
elems(2).SampleTime = -1;
elems(2).Complexity = 'real';
elems(2).Min = [];
elems(2).Max = [];
elems(2).DocUnits = '';
elems(2).Description = '';

deltaeBus = Simulink.Bus;
deltaeBus.HeaderFile = '';
deltaeBus.Description = '';
deltaeBus.DataScope = 'Auto';
deltaeBus.Alignment = -1;
deltaeBus.Elements = elems;
clear elems;
assignin('base','deltaeBus', deltaeBus);

% Bus object: deltarBus 
clear elems;
elems(1) = Simulink.BusElement;
elems(1).Name = 'idxdeltar';
elems(1).Dimensions = 1;
elems(1).DimensionsMode = 'Fixed';
elems(1).DataType = 'int32';
elems(1).SampleTime = -1;
elems(1).Complexity = 'real';
elems(1).Min = [];
elems(1).Max = [];
elems(1).DocUnits = '';
elems(1).Description = '';

elems(2) = Simulink.BusElement;
elems(2).Name = 'fdeltar';
elems(2).Dimensions = 1;
elems(2).DimensionsMode = 'Fixed';
elems(2).DataType = 'double';
elems(2).SampleTime = -1;
elems(2).Complexity = 'real';
elems(2).Min = [];
elems(2).Max = [];
elems(2).DocUnits = '';
elems(2).Description = '';

deltarBus = Simulink.Bus;
deltarBus.HeaderFile = '';
deltarBus.Description = '';
deltarBus.DataScope = 'Auto';
deltarBus.Alignment = -1;
deltarBus.Elements = elems;
clear elems;
assignin('base','deltarBus', deltarBus);

% Bus object: deltafrBus 
clear elems;
elems(1) = Simulink.BusElement;
elems(1).Name = 'idxdeltafr';
elems(1).Dimensions = 1;
elems(1).DimensionsMode = 'Fixed';
elems(1).DataType = 'int32';
elems(1).SampleTime = -1;
elems(1).Complexity = 'real';
elems(1).Min = [];
elems(1).Max = [];
elems(1).DocUnits = '';
elems(1).Description = '';

elems(2) = Simulink.BusElement;
elems(2).Name = 'fdeltafr';
elems(2).Dimensions = 1;
elems(2).DimensionsMode = 'Fixed';
elems(2).DataType = 'double';
elems(2).SampleTime = -1;
elems(2).Complexity = 'real';
elems(2).Min = [];
elems(2).Max = [];
elems(2).DocUnits = '';
elems(2).Description = '';

deltafrBus = Simulink.Bus;
deltafrBus.HeaderFile = '';
deltafrBus.Description = '';
deltafrBus.DataScope = 'Auto';
deltafrBus.Alignment = -1;
deltafrBus.Elements = elems;
clear elems;
assignin('base','deltafrBus', deltafrBus);

% Bus object: deltaflBus 
clear elems;
elems(1) = Simulink.BusElement;
elems(1).Name = 'idxdeltafl';
elems(1).Dimensions = 1;
elems(1).DimensionsMode = 'Fixed';
elems(1).DataType = 'int32';
elems(1).SampleTime = -1;
elems(1).Complexity = 'real';
elems(1).Min = [];
elems(1).Max = [];
elems(1).DocUnits = '';
elems(1).Description = '';

elems(2) = Simulink.BusElement;
elems(2).Name = 'fdeltafl';
elems(2).Dimensions = 1;
elems(2).DimensionsMode = 'Fixed';
elems(2).DataType = 'double';
elems(2).SampleTime = -1;
elems(2).Complexity = 'real';
elems(2).Min = [];
elems(2).Max = [];
elems(2).DocUnits = '';
elems(2).Description = '';

deltaflBus = Simulink.Bus;
deltaflBus.HeaderFile = '';
deltaflBus.Description = '';
deltaflBus.DataScope = 'Auto';
deltaflBus.Alignment = -1;
deltaflBus.Elements = elems;
clear elems;
assignin('base','deltaflBus', deltaflBus);

