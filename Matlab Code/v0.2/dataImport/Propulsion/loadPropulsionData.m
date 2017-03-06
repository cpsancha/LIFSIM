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


%% MAP THRUST, POWER, CURRENT AND RPM IN FUNCTION OF VOLTAGE, HEIGHT AND VELOCITY
% Mapping folder in which are stored the generated files
filename = mfilename('fullpath');
[pathstr,~,~] = fileparts(filename);
mappingFolder = 'motorPropMapping';
mappingPath = strcat(pathstr,filesep,mappingFolder);

% Ranges:
voltageLength =  60;
heightLength  =  10;
vflightLength =  90;
    
% Limits:
voltageLim = [0,   30]; %V
heightLim  = [0, 1000]; %m
vflightLim = [0,   85]; %m/s
    
% Mapping:    
if isequal(length(motorNames),length(propellerNames))
    
    %Store files in the mapping directory as list
    fileList = dir(mappingPath);
    
    for i=1:length(motorNames)
        %Create expected file filename
        motorModel     = LD.Propulsion.(motorNames{i}).Model;
        propellerModel = LD.Propulsion.(propellerNames{i}).Model;
        expectedFile   = strcat(motorModel,propellerModel);
        expectedFile = regexprep(expectedFile,'[\\\*:/><|?~!@#$%^&()_\-{}\s]','');
        expectedSufix = '.mat';
        index = zeros(1,size(fileList,1));
        
        %Check if expected file exists
        for j=1:size(fileList,1)
            index(j) = isequal(char(fileList(j).name),strcat(expectedFile,expectedSufix));
        end
        
        if sum(index)
            %Read mapped data from file
            [validData,Voltage,Height,Vflight,Thrust,Power,Current,RPM,exitFlag] = ...
                loadMappedData(strcat(expectedFile,expectedSufix),motorModel,...
                propellerModel,voltageLim,heightLim,vflightLim,voltageLength,...
                heightLength,vflightLength);
            
            if validData
                %Store mapped data in motorData struct
                LD.Propulsion.(motorNames{i}).Voltage  = Voltage;
                LD.Propulsion.(motorNames{i}).Height   = Height;
                LD.Propulsion.(motorNames{i}).Vflight  = Vflight;
                LD.Propulsion.(motorNames{i}).Thrust   = Thrust;
                LD.Propulsion.(motorNames{i}).Power    = Power;
                LD.Propulsion.(motorNames{i}).Current  = Current;
                LD.Propulsion.(motorNames{i}).RPM      = RPM;
                LD.Propulsion.(motorNames{i}).exitFlag = exitFlag;
            else
                %Calculate mapped data
                    [Voltage,Height,Vflight,Thrust,Power,Current,RPM,exitFlag] = ...
                     propulsionDataMapping(LD.Propulsion.(motorNames{i}),LD...
                     .Propulsion.(propellerNames{i}),voltageLim,heightLim,...
                     vflightLim,voltageLength,heightLength,vflightLength);
            
                %Save mapped data in mat file
                save(strcat(mappingPath,filesep,expectedFile,expectedSufix),...
                     'motorModel','propellerModel',...
                     'voltageLim','heightLim','vflightLim',...
                     'voltageLength','heightLength','vflightLength',...
                     'Thrust','Power','Current','RPM','exitFlag')
                
                %Store mapped data in motorData struct
                LD.Propulsion.(motorNames{i}).Voltage  = Voltage;
                LD.Propulsion.(motorNames{i}).Height   = Height;
                LD.Propulsion.(motorNames{i}).Vflight  = Vflight;
                LD.Propulsion.(motorNames{i}).Thrust   = Thrust;
                LD.Propulsion.(motorNames{i}).Power    = Power;
                LD.Propulsion.(motorNames{i}).Current  = Current;
                LD.Propulsion.(motorNames{i}).RPM      = RPM;
                LD.Propulsion.(motorNames{i}).exitFlag = exitFlag;
            end
            
        else
            %Calculate mapped data
            [Voltage,Height,Vflight,Thrust,Power,Current,RPM,exitFlag] = ...
                 propulsionDataMapping(LD.Propulsion.(motorNames{i}),LD...
                 .Propulsion.(propellerNames{i}),voltageLim,heightLim,...
                 vflightLim,voltageLength,heightLength,vflightLength);
            
            %Save mapped data in mat file
            save(strcat(mappingPath,filesep,expectedFile,expectedSufix),...
                 'motorModel','propellerModel',...
                 'voltageLim','heightLim','vflightLim',...
                 'voltageLength','heightLength','vflightLength',...
                 'Thrust','Power','Current','RPM','exitFlag')
            
            %Store mapped data in motorData struct
            LD.Propulsion.(motorNames{i}).Voltage  = Voltage;
            LD.Propulsion.(motorNames{i}).Height   = Height;
            LD.Propulsion.(motorNames{i}).Vflight  = Vflight;
            LD.Propulsion.(motorNames{i}).Thrust   = Thrust;
            LD.Propulsion.(motorNames{i}).Power    = Power;
            LD.Propulsion.(motorNames{i}).Current  = Current;
            LD.Propulsion.(motorNames{i}).RPM      = RPM;
            LD.Propulsion.(motorNames{i}).exitFlag = exitFlag;
            
        end
    end
else
    wrn = msgbox('The number of motors is not consistent with the number of propellers.','Warning','warn');
        uiwait(wrn);
        disp('Inconsistent motor-propeller dimensions, please fix me...')
        pause
end




%% DEFINE ESC
LD.Propulsion.ESC1.eta = 0.99;
LD.Propulsion.ESC2.eta = 0.99;
LD.Propulsion.ESC3.eta = 0.99;
LD.Propulsion.ESC4.eta = 0.99;
LD.Propulsion.ESC5.eta = 0.99;



%% DEFINE INITIAL THROTTLE
LD.Propulsion.Throttle1 = 1;
LD.Propulsion.Throttle2 = 1;
LD.Propulsion.Throttle3 = 1;
LD.Propulsion.Throttle4 = 1;
LD.Propulsion.Throttle5 = 1;




%% CLEAR USED VARIABLES
clear motorData motorNames motorModels index wrn i j propellerData motorPosition
clear propellerData propellerFiles propellerNames propellerModels
clear fileList filename pathstr mappingFolder mappingPath expectedFile expectedSufix
clear Voltage Height Vflight voltageLim voltageLength heightLim heightLength vflightLim vflightLength
clear motorModel propellerModel validData Thrust Power Current RPM exitFlag    