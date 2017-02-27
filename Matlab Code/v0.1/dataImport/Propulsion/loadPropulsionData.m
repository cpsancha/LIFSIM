%% VARIABLE DEFINITIONS:    
%          J=V/nD (advance ratio)                                                                               
%          Ct=T/(rho * n**2 * D**4) (thrust coef.)                                                              
%          Cp=P/(rho * n**3 * D**5) (power coef.)                                                               
%          Pe=Ct*J/Cp (efficiency)                                                                              
%          V          (model speed in MPH)             
%          PWR        (Hp)   
%          Torque     (In-Lbf)     
%          Thrust     (Lbf) 


%% DEFINE MOTORS

%Import motor data library
motorData = importMotorData();
    
%Define motors' names in the structure
% Motor order is the following: Horizontal flight, Front right wing, 
% Rear right wing, Rear left wing, Front left wing
motorNames = string({'motor1',...
                     'motor2',...
                     'motor3',...
                     'motor4',...
                     'motor5'});

%Define a matrix of motors' propeller position in the same order as above
% Each row is the position of the motor's propeller [x y z] from the nose in meters
motorPosition = [[-2.012  0.000 0.000];...
                 [-0.643  1.151 0.000];...
                 [-1.440  1.151 0.000];...
                 [-1.440 -1.151 0.000];...
                 [-0.643 -1.151 0.000]];

%Define motors' models in the same order as above    
motorModels = string({'AXI 5320/28 GOLD LINE V2',...
                      'AXI 5320/28 GOLD LINE V2',...
                      'AXI 5320/28 GOLD LINE V2',...
                      'AXI 5320/28 GOLD LINE V2',...
                      'AXI 5320/28 GOLD LINE V2'});
    
%Define motors in LD structure
for i=1:length(motorNames)
    if sum(contains(string({motorData.Model}), motorModels{i}))
        index = find(contains(string({motorData.Model}),motorModels{i}));
        LD.Propulsion.(motorNames{i}) = motorData(index); %#ok<FNDSB>
        LD.Propulsion.(motorNames{i}).Position = motorPosition(i,:); 
    else
        eval(sprintf('wrn = msgbox(''Motor selection for %s is not present in the database.'',''Warning'',''warn'');',motorNames{i}));
        uiwait(wrn);
        disp('Wrong motor selection, please fix me...')
        pause
    end
end


%% DEFINE PROPELLERS 

%Define the propellers' files from APC to be imported into the database
% More files can be found in the following web page:
% https://www.apcprop.com/v/PERFILES_WEB/listDatafiles.asp
propellerFiles = string({'PER3_17x8E.txt',...
                         'PER3_15x10E.txt'});

%Import propeller data library
for i=1:length(propellerFiles)
    propellerData(i) = importPropellerData(propellerFiles{i}); %#ok<SAGROW>
end

%Define propellers' names in the structure
% Propeller order is the following: Horizontal flight, Front right wing, 
% Rear right wing, Rear left wing, Front left wing
propellerNames = string({'propeller1',...
                         'propeller2',...
                         'propeller3',...
                         'propeller4',...
                         'propeller5'});
                     
%Define propellers' models in the same order as above    
propellerModels = string({'15x10E',...
                          '17x8E',...
                          '17x8E',...
                          '17x8E',...
                          '17x8E'});
                  
%Define propellers in LD structure
for i=1:length(propellerNames)
    if sum(contains(string({propellerData.Model}), propellerModels{i}))
        index = find(contains(string({propellerData.Model}),propellerModels{i}));
        LD.Propulsion.(propellerNames{i}) = propellerData(index); %#ok<FNDSB>
    else
        eval(sprintf('wrn = msgbox(''Propeller selection for %s is not present in the database.'',''Warning'',''warn'');',propellerNames{i}));
        uiwait(wrn);
        disp('Wrong propeller selection, please fix me...')
        pause
    end
end                  




%% DEFINE INITIAL THROTTLE
LD.Propulsion.Throttle1 = 0.5;
LD.Propulsion.Throttle2 = 0.5;
LD.Propulsion.Throttle3 = 0.5;
LD.Propulsion.Throttle4 = 0.5;
LD.Propulsion.Throttle5 = 0.5;


%% CLEAR USED VARIABLES
clear motorData motorNames motorModels index wrn i propellerData motorPosition
clear propellerData propellerFiles propellerNames propellerModels
    