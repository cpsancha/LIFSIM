%% IMPORT LIBIS DATA TO DATCOM 1976 STRUCTURE
%Script to import the static and dynamic data of the LIBIS RPAS in the
%DATCOM structure style 1976 Version (File Type 6)
%     load('LIBISData.mat')
    LD.case = 'LIBIS RPAS';

    
%% DEFINE UNITS
%Define the used units in the inputs
% CHECK. Validar si escoger correctamente las unidades de entrada
% afecta al módulo y que nuestras derivadas estén en m, rad
    LD.dim   = 'M';   %Character vector denoting the specified system 
                      %of units for the case (FT(default), IN, M or CM)
    LD.deriv = 'RAD'; %Character vector denoting the specified angle 
                      %units for the case (DEG(default) or RAD)


%% SET THE REFERENCE DIMENSIONS
%Dimensiones de referencia (excel page: "Derivadas estabilidad Long.)
    LD.sref  = 1.251;  %Scalar denoting the reference area (2 wings)
    LD.cbar  = 0.284;  %Scalar denoting the longitudinal reference length (1 wing)
    LD.blref = 2.237;  %Scalar denoting the lateral reference length (1 wing)
    LD.xcg   = 1.015;  %Distance from the nose of the plane to center of gravity (meters)
    LD.xba   = 0.744;  %Distance from the nose of the plane to the wing's leading edge (meters)
    LD.xcg_CMA = (LD.xcg - LD.xba)/LD.cbar; %Position of the xcg in the aerodynamic mean chord

    
    
%% SET THE CASES OF ANALYSIS
%Stability derivatives vary with the following parameters. Parse parameters
%with values for which stability derivatives information is avaliable.
%NEVER --> Leave empty array if no variation is considered for the parameter.
    %Array of angles of attack.
%         LD.alpha  = [3.03*pi/180,9];
        LD.alpha  = [1,2];
    %Array of altitudes.
        LD.alt    = [91.44,100];
    %Array of positions of the center of gravity defined as fraction of the CMA (0<xcg<1). 
        LD.xcg  = [0.955,1];
    %Array of left and right lifting surface streamwise control deflection 
    %angles, which are defined positive for trailing-edge down.
        LD.deltae  = [0,1];
        LD.deltar  = [0,1];
        LD.deltafr = [0,1];
        LD.deltafl = [0,1];
    %Array of Mach numbers
        LD.mach   = [0.12];
    %Array of Rynolds numbers
        LD.rnnub  = [];
    %Array of control-surface streamwise deflection angles
        LD.delta  = [];

        
        
  %% READ THE LENGTH OF THE VECTORS
  %Store the length of the cases of analysis
        LD.nmach  = length(LD.mach);
        LD.nalt   = length(LD.alt);
        LD.nalpha = length(LD.alpha);
        LD.nxcg   = length(LD.xcg);
        %LD.ndelta=length(LD.delta); Está mal, sería de 0 a 3 en función de
        %si existen delta, deltal y  deltar
        
        
%% 
%Stability derivatives matrix initialization:
%Static
% LD.cd=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build),length(LD.grndht),length(LD.delta));%
% LD.cl=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build),length(LD.grndht),length(LD.delta));%
% LD.cm=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build),length(LD.grndht),length(LD.delta));%
% LD.cn=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build),length(LD.grndht),length(LD.delta));
% LD.ca=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build),length(LD.grndht),length(LD.delta));
% LD.xcp=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build),length(LD.grndht),length(LD.delta));
% LD.cla=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build),length(LD.grndht),length(LD.delta));
% LD.cma=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build),length(LD.grndht),length(LD.delta));
% LD.cyb=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build),length(LD.grndht),length(LD.delta));%
% LD.cnb=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build),length(LD.grndht),length(LD.delta));%
% LD.clb=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build),length(LD.grndht),length(LD.delta));%
% %Dynamic
% LD.clq=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build));%
% LD.cmq=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build));%
% LD.clad=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build));%
% LD.cmad=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build));%
% LD.clp=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build));%
% LD.cyp=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build));%
% LD.cnp=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build));%
% LD.cnr=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build));%
% LD.clr=zeros(length(LD.alpha),length(LD.mach),length(LD.alt),...
%     length(LD.build));%
% %High Lift (Not applicable but initialized)
% LD.dcl_sym=zeros(length(LD.alpha),length(LD.mach),length(LD.alt));
% LD.dcm_sym=zeros(length(LD.alpha),length(LD.mach),length(LD.alt));
% LD.dclmax_sym=zeros(length(LD.alpha),length(LD.mach),length(LD.alt));
% LD.dcdmin_sym=zeros(length(LD.alpha),length(LD.mach),length(LD.alt));
% LD.clad_sym=zeros(length(LD.alpha),length(LD.mach),length(LD.alt));
% LD.cha_sym=zeros(length(LD.alpha),length(LD.mach),length(LD.alt));
% LD.chd_sym=zeros(length(LD.alpha),length(LD.mach),length(LD.alt));
% LD.dcid_sym=zeros(length(LD.alpha),length(LD.mach),length(LD.alt));
%         

        
        
        
        
        


