%Script to import the static and dynamic data of the LIBIS RPAS in the
%DATCOM structure style 1976 Version (File Type 6)

% load('LIBISData.mat')


%Stability derivatives vary with the following parameters.Parse parameters
%with values for which stability derivatives information is avaliable.
%Leave empty array if no variation is considered for the parameter.

    %Array of angles of attack.
        LD.alpha  = [3.03*pi/180];
    %Array of Mach numbers
        LD.mach   = [0.12];
    %Array of altitudes.
        LD.alt    = [91.44];
    %Array of configurations. Will be used for changes in Xcg
        LD.build  = [];
    %Array of ground heights.
        LD.grndht = [];
    %Array of control-surface streamwise deflection angles
        LD.delta  = [];  
        
        
        LD.nmach=length(LD.mach);
        LD.nalt=length(LD.alt);
        LD.alpha=length(LD.alpha);
        LD.ndelta=length(LD.delta);
        LD.ngh=length(LD.grndht);
        
        %% Dimensiones de referencia (excel page: "Derivadas estabilidad Long.)
        
        LD.sref=0.63;  %reference area in m  (2 wings)
        LD.cbar=0.28;  %longitudinal reference length (1 wing)
        LD.blref=2.24; %lateral reference length (1 wing)
        
        %% CHECK. validar si escoger correctamente las unidades de entrada
        %afecta al módulo y que nuestras derivadas estén en m, rad
        
        LD.dim='m'; %Character vector denoting the specified system of units for the case
        LD.deriv='rad'; %Character vector denoting the specified angle units for the case
        
        
        


