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
    LD.xba   = -0.744;  %Distance from the nose of the plane to the wing's leading edge (meters),
    %in body axis with origin at the nose.
%   LD.xcg_CMA = (LD.xcg - LD.xba)/LD.cbar; %Position of the xcg in the aerodynamic mean chord

    
%% CASES OF ANALYSIS AND COEFFICIENTS DEFINITION
LD.Stability.analysisCases = string({'alpha','beta','alt','xcg','deltae','deltar','deltafr','deltafl'});
LD.Stability.coeffs = string({'CD0','CDalpha','CDalpha_dot','CDq','CDdeltae','CDdeltafr','CDdeltafl',...
                        'CL0','CLalpha','CLalpha_dot','CLq','CLdeltae','CLdeltafr','CLdeltafl',...
                        'Cm0','Cmalpha','Cmalpha_dot','Cmq','Cmdeltae','Cmdeltafr','Cmdeltafl',...
                        'CYbeta','CYp','CYr','CYdeltar','CYdeltafr','CYdeltafl',...
                        'Clbeta','Clp','Clr','Cldeltar','Cldeltafr','Cldeltafl',...
                        'Cnbeta','Cnp','Cnr','Cndeltar','Cndeltafr','Cndeltafl'});
                    
                    
%% SET THE CASES OF ANALYSIS
%Stability derivatives vary with the following parameters. Parse parameters
%with values for which stability derivatives information is avaliable.
%NEVER --> Leave empty array if no variation is considered for the parameter.
    %Array of angles of attack in radians.
        LD.alpha = [-180, -90, -12.74, -12, 0, 12, 12.74, 90, 180]; %In degrees for convenience
        LD.alpha = LD.alpha.*pi/180; %In radians
    %Array of sideslip angles
        LD.beta = 0; %In degrees for convenience
        LD.beta = LD.beta.*pi/180; %In radians
    %Array of altitudes.
        LD.alt    = 91.44;
    %Array of positions of the center of gravity. Distance from the nose of
    %the plane to center of gravity (meters) in body axis with origin at the
    %nose
        LD.xcg    = -1.013;
    %Array of elevator, rudder, and right and left flaperons streamwise control 
    %deflection angles, which are defined positive for trailing-edge down.
        LD.deltae  = 10*pi/180;
        LD.deltar  = 10*pi/180;
        LD.deltafr = 10*pi/180;
        LD.deltafl = 10*pi/180;



        
        
  %% READ THE LENGTH OF THE VECTORS
  %Store the length of the cases of analysis
%         LD.nalpha = length(LD.alpha);
%         LD.nalt   = length(LD.alt);
%         LD.nxcg   = length(LD.xcg);

        
        
%% LONGITUDINAL STABILITY DERIVATIVES
%Parasite coeffs
    %CD0
    CD0 = [0, 0, 0, 0.024, 0.024, 0.024, 0, 0, 0];
    LD.Stability.CD0 = fillStabilityCoeffs(LD, string({'alpha'}), CD0);
    %CL0
    CL0 = [0, 0, 0, 0.181, 0.181, 0.181, 0, 0, 0];
    LD.Stability.CL0 = fillStabilityCoeffs(LD, string({'alpha'}), CL0);
    %Cm0
    Cm0 = [0, 0, 0, -0.003, -0.003, -0.003, 0, 0, 0];
    LD.Stability.Cm0 = fillStabilityCoeffs(LD, string({'alpha'}), Cm0);
    
%Alpha
    %CDalpha
    CDalpha = [0, 0, 0, 0.239, 0.239, 0.239, 0, 0, 0];
    LD.Stability.CDalpha = fillStabilityCoeffs(LD, string({'alpha'}), CDalpha);
    %CLalpha
    CLalpha = [0, 0, 0, 4.138, 4.138, 4.138, 0, 0, 0];
    LD.Stability.CLalpha = fillStabilityCoeffs(LD, string({'alpha'}), CLalpha);
    %Cmalpha
    Cmalpha = [0, 0, 0, -0.7074, -0.707, -0.707, 0, 0, 0];
    LD.Stability.Cmalpha = fillStabilityCoeffs(LD, string({'alpha'}), Cmalpha);
    
%Alpha_dot
    %CDalpha_dot
    CDalpha_dot = [0, 0, 0, 0.000, 0.000, 0.000, 0, 0, 0];
    LD.Stability.CDalpha_dot = fillStabilityCoeffs(LD, string({'alpha'}), CDalpha_dot);
    %CLalpha_dot
    CLalpha_dot = [0, 0, 0, 7.411, 7.411, 7.411, 0, 0, 0];
    LD.Stability.CLalpha_dot = fillStabilityCoeffs(LD, string({'alpha'}), CLalpha_dot);
    %Cmalpha_dot
    Cmalpha_dot = [0, 0, 0, -4.063, -4.063, -4.063, 0, 0, 0];
    LD.Stability.Cmalpha_dot = fillStabilityCoeffs(LD, string({'alpha'}), Cmalpha_dot);
    
%q
    %CDq
    CDq = [0, 0, 0, 0.000, 0.000, 0.000, 0, 0, 0];
    LD.Stability.CDq = fillStabilityCoeffs(LD, string({'alpha'}), CDq);
    %CLq
    CLq = [0, 0, 0, 6.248, 6.248, 6.248, 0, 0, 0];
    LD.Stability.CLq = fillStabilityCoeffs(LD, string({'alpha'}), CLq);
    %Cmq
    Cmq = [0, 0, 0, -18.410, -18.410, -18.410, 0, 0, 0];
    LD.Stability.Cmq = fillStabilityCoeffs(LD, string({'alpha'}), Cmq);
    
%deltae --> elevator
    %CDdeltae
    CDdeltae = [0, 0, 0, 0.005, 0.005, 0.005, 0, 0, 0];
    LD.Stability.CDdeltae = fillStabilityCoeffs(LD, string({'alpha'}), CDdeltae);
    %CLdeltae
    CLdeltae = [0, 0, 0, 1.391, 1.391, 1.391, 0, 0, 0];
    LD.Stability.CLdeltae = fillStabilityCoeffs(LD, string({'alpha'}), CLdeltae);
    %Cmdeltae
    Cmdeltae = [0, 0, 0, -2.930, -2.930, -2.930, 0, 0, 0];
    LD.Stability.Cmdeltae = fillStabilityCoeffs(LD, string({'alpha'}), Cmdeltae);
    
%deltafr --> rigth flaperon
    %CDdeltafr
    CDdeltafr = LD.Stability.CDdeltae/2; %TO BE CALCULATED
    LD.Stability.CDdeltafr = fillStabilityCoeffs(LD, string({'alpha'}), CDdeltafr);
    %CLdeltafr
    CLdeltafr = LD.Stability.CLdeltae/2; %TO BE CALCULATED
    LD.Stability.CLdeltafr = fillStabilityCoeffs(LD, string({'alpha'}), CLdeltafr);
    %Cmdeltafr
    Cmdeltafr = [0, 0, 0, 0.000, 0.000, 0.000, 0, 0, 0]; %TO BE CALCULATED
    LD.Stability.Cmdeltafr = fillStabilityCoeffs(LD, string({'alpha'}), Cmdeltafr);
    
%deltafl --> left flaperon
    %CDdeltafl
    CDdeltafl = LD.Stability.CDdeltae/2; %TO BE CALCULATED
    LD.Stability.CDdeltafl = fillStabilityCoeffs(LD, string({'alpha'}), CDdeltafl);
    %CLdeltafl
    CLdeltafl = LD.Stability.CLdeltae/2; %TO BE CALCULATED
    LD.Stability.CLdeltafl = fillStabilityCoeffs(LD, string({'alpha'}), CLdeltafl);  
    %Cmdeltafl
    Cmdeltafl = [0, 0, 0, 0.000, 0.000, 0.000, 0, 0, 0]; %TO BE CALCULATED
    LD.Stability.Cmdeltafl = fillStabilityCoeffs(LD, string({'alpha'}), Cmdeltafl);
    

%% LATERAL-DIRECTIONAL DERIVATIVES
%beta
    %CYbeta
    CYbeta = [0, 0, 0, -0.4250, -0.4250, -0.4250, 0, 0, 0];
    LD.Stability.CYbeta = fillStabilityCoeffs(LD, string({'alpha'}), CYbeta);
    %Clbeta
    Clbeta = [0, 0, 0, -0.0286, -0.0286, -0.0286, 0, 0, 0];
    LD.Stability.Clbeta = fillStabilityCoeffs(LD, string({'alpha'}), Clbeta);
    %Cnbeta
    Cnbeta = [0, 0, 0, 0.0170, 0.0170, 0.0170, 0, 0, 0];
    LD.Stability.Cnbeta = fillStabilityCoeffs(LD, string({'alpha'}), Cnbeta);
    
%p
    %CYp
    CYp = [0, 0, 0, -0.004, -0.004, -0.004, 0, 0, 0];
    LD.Stability.CYp = fillStabilityCoeffs(LD, string({'alpha'}), CYp);
    %Clp
    Clp = [0, 0, 0, -0.472, -0.472, -0.472, 0, 0, 0];
    LD.Stability.Clp = fillStabilityCoeffs(LD, string({'alpha'}), Clp);
    %Cnp
    Cnp = [0, 0, 0, -0.019, -0.019, -0.019, 0, 0, 0];
    LD.Stability.Cnp = fillStabilityCoeffs(LD, string({'alpha'}), Cnp);
    
%r
    %CYr
    CYr = [0, 0, 0, 0.120, 0.120, 0.120, 0, 0, 0];
    LD.Stability.CYr = fillStabilityCoeffs(LD, string({'alpha'}), CYr);
    %Clr
    Clr = [0, 0, 0, 0.049, 0.049, 0.049, 0, 0, 0];
    LD.Stability.Clr = fillStabilityCoeffs(LD, string({'alpha'}), Clr);
    %Cnr
    Cnr = [0, 0, 0, -0.040, -0.040, -0.040, 0, 0, 0];
    LD.Stability.Cnr = fillStabilityCoeffs(LD, string({'alpha'}), Cnr);
    
%deltar --> rudder
    %CYdeltar
    CYdeltar = [0, 0, 0, -0.212, -0.212, -0.212, 0, 0, 0];
    LD.Stability.CYdeltar = fillStabilityCoeffs(LD, string({'alpha'}), CYdeltar);
    %Cldeltar
    Cldeltar = [0, 0, 0, -0.002, -0.002, -0.002, 0, 0, 0];
    LD.Stability.Cldeltar = fillStabilityCoeffs(LD, string({'alpha'}), Cldeltar);
    %Cndeltar
    Cndeltar = [0, 0, 0, 0.051, 0.051, 0.051, 0, 0, 0];
    LD.Stability.Cndeltar = fillStabilityCoeffs(LD, string({'alpha'}), Cndeltar);
    
%deltafr --> rigth flaperon
    %CYdeltafr
    CYdeltafr = [0, 0, 0, 0.000, 0.000, 0.000, 0, 0, 0];
    LD.Stability.CYdeltafr = fillStabilityCoeffs(LD, string({'alpha'}), CYdeltafr);
    %Cldeltafr
    Cldeltafr = [0, 0, 0, -0.125, -0.125, -0.125, 0, 0, 0]; %-deltaa/2
    LD.Stability.Cldeltafr = fillStabilityCoeffs(LD, string({'alpha'}), Cldeltafr);
    %Cndeltafr
    Cndeltafr = [0, 0, 0, 0.001, 0.001, 0.001, 0, 0, 0]; %-deltaa/2
    LD.Stability.Cndeltafr = fillStabilityCoeffs(LD, string({'alpha'}), Cndeltafr); 
    
%deltafl --> left flaperon
    %CYdeltafl
    CYdeltafl = [0, 0, 0, 0.000, 0.000, 0.000, 0, 0, 0];
    LD.Stability.CYdeltafl = fillStabilityCoeffs(LD, string({'alpha'}), CYdeltafl);
    %Cldeltafl
    Cldeltafl = [0, 0, 0, 0.125, 0.125, 0.125, 0, 0, 0]; %deltaa/2
    LD.Stability.Cldeltafl = fillStabilityCoeffs(LD, string({'alpha'}), Cldeltafl); 
    %Cndeltafl
    Cndeltafl = [0, 0, 0, -0.001, -0.001, -0.001, 0, 0, 0]; %deltaa/2
    LD.Stability.Cndeltafl = fillStabilityCoeffs(LD, string({'alpha'}), Cndeltafl); 
        
        
        
%% CLEAR USED VARIABLES
clear CD0 CDalpha CDalpha_dot CDq CDdeltae CDdeltafr CDdeltafl
clear CL0 CLalpha CLalpha_dot CLq CLdeltae CLdeltafr CLdeltafl
clear Cm0 Cmalpha Cmalpha_dot Cmq Cmdeltae Cmdeltafr Cmdeltafl
clear CYbeta CYp CYr CYdeltar CYdeltafr CYdeltafl
clear Clbeta Clp Clr Cldeltar Cldeltafr Cldeltafl
clear Cnbeta Cnp Cnr Cndeltar Cndeltafr Cndeltafl                  

