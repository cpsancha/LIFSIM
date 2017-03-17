function tests = linearTest()
% linearTest - Sample test file with one test point
% "testLinearComparison" that makes sure that the outputs of the model for
% ss the linear model (linsys) of the aircraft and nonlinear model are
% equal or within a specific tolerance. To run the test simply run
% "runtests('linearTest')" from the command window or runProjectTests.

% Copyright 2013 The MathWorks, Inc.

tests = functiontests(localfunctions);

end

function testLinearComparison(testCase) %#okDEFNU
    %% Load Model
    model = 'trimNonlinearAirframe';
    load_system(model);

    %% Set Variants
    oldVariants = evalin('base','Variants');
    Variants = oldVariants;
    Variants.Vehicle = 1; %Set Vehicle to Nonlinear
    Variants.Actuators = 3;
    assignin('base','Variants',Variants);

    %% Load linear model trim point
    load('linearizedAirframe.mat');
    assignin('base','op',op);
    
    %% Load maneuver for testing and map it to the inputs
    load('maneuverHover.mat');
    assignin('base','maneuverHover',maneuverHover);
    inputMap = getSlRootInportMap('model',model,'MappingMode','SignalName',...
    'SignalValue',{maneuverHover},'SignalName',{'maneuverHover'});
    inputString = getInputString(inputMap,'base');

    %% Set initial states and inputs
    set_param(model,'LoadExternalInput','on');
    set_param(model,'ExternalInput',inputString);
    set_param(model,'LoadInitialState','on');
    set_param(model,'InitialState','getstatestruct(op)');
    set_param(model,'SaveOutput','on');

    %% Simulate Nonlinear model using trim conditions
    sim(model);

    %% Assign expected values for final time step
    expSimVb = yout.signals(1).values;
    expSimOmega = yout.signals(2).values;

    %% Change variant to linear
    Variants.Vehicle = 0;
    assignin('base','Variants',Variants);
    set_param(model,'LoadInitialState','off');    
    
    %% Simulate Linear model using trim conditions and maneuver
    sim(model);

    %% Assign actual values for final time step
    actSimVb = yout.signals(1).values;
    actSimOmega = yout.signals(2).values;

    %% Run test point
    verifyEqual(testCase,actSimVb,expSimVb,'AbsTol',0.05);
    verifyEqual(testCase,actSimOmega,expSimOmega,'AbsTol',1e-10);

    %% Clean up
    assignin('base','Variants',oldVariants);
    close_system(model,0);
    bdclose({'linearAirframe','nonlinearAirframe'});
    evalin('base', 'clear(''op'',''opreport'',''linsys'')');
end