%% IMPORT LIBIS DATA TO DATCOM 1976 STRUCTURE
%Script to import the static and dynamic data of the LIBIS RPAS in the
%DATCOM structure style 1976 Version (File Type 6)
    % load('LIBISData.mat')
    LD.case = 'LIBIS RPAS';
    LD.version = 1976;

    
%% DEFINE UNITS
%Define the used units in the inputs
% CHECK. Validar si escoger correctamente las unidades de entrada
% afecta al módulo y que nuestras derivadas estén en m, rad
    LD.dim   = 'M';   %Character vector denoting the specified system 
                      %of units for the case (FT(default), IN, M or CM)
    LD.deriv = 'RAD'; %Character vector denoting the specified angle 
                      %units for the case (DEG(default) or RAD)

    
%% DEFINE OPTIONS
%Logical denoting the reading of dynamic derivative data for the case.     
%When dynamic derivative runs are read, this value is set to true.
	LD.damp = true; %False by default


%% SET THE REFERENCE DIMENSIONS
%Dimensiones de referencia (excel page: "Derivadas estabilidad Long.)
    LD.sref  = 0.63;  %Scalar denoting the reference area (2 wings)
    LD.cbar  = 0.28;  %Scalar denoting the longitudinal reference length (1 wing)
    LD.blref = 2.24;  %Scalar denoting the lateral reference length (1 wing)
    
    
%% SET THE CASES OF ANALYSIS
%Stability derivatives vary with the following parameters. Parse parameters
%with values for which stability derivatives information is avaliable.
%Leave empty array if no variation is considered for the parameter.
    %Array of angles of attack.
        LD.alpha  = [3.03*pi/180];  %¿¿RADIANS???
    %Array of Mach numbers
        LD.mach   = [0.12];
    %Array of altitudes.
        LD.alt    = [91.44];
    %Array of Rynolds numbers
        LD.rnnub  = [];
    %Array of configurations. Will be used for changes in Xcg 
%FIXME No se usa así
        LD.build  = [];
    %Array of ground heights.
        LD.grndht = [];
    %Array of control-surface streamwise deflection angles
        LD.delta  = [];
    %Array of left and right lifting surface streamwise control deflection 
    %angles, which are defined positive for trailing-edge down.
        LD.deltal  = [];
        LD.deltar  = [];
        
        
  %% READ THE LENGTH OF THE VECTORS
  %Store the length of the cases of analysis
        LD.nmach=length(LD.mach);
        LD.nalt=length(LD.alt);
        LD.nalpha=length(LD.alpha);
        %LD.ndelta=length(LD.delta); Está mal, sería de 0 a 3 en función de
        %si existen delta, deltal y  deltar
        LD.ngh=length(LD.grndht);
        

        
        
        
        
        


