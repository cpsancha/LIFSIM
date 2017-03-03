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

    
    
%% SET THE CASES OF ANALYSIS
%Stability derivatives vary with the following parameters. Parse parameters
%with values for which stability derivatives information is avaliable.
%NEVER --> Leave empty array if no variation is considered for the parameter.
    %Array of angles of attack in radians.
        LD.alpha  = 3.03*pi/180;
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
    LD.Stability.CD0 =  0.024;
    LD.Stability.CL0 =  0.181;
    LD.Stability.Cm0 = -0.003;
%Alpha
    LD.Stability.CDalpha =  0.239;
    LD.Stability.CLalpha =  4.138;
    LD.Stability.Cmalpha = -0.707;
%Alpha_dot
    LD.Stability.CDalpha_dot =  0.000;
    LD.Stability.CLalpha_dot =  7.411;
    LD.Stability.Cmalpha_dot = -4.063;
%q
    LD.Stability.CDq =   0.000;
    LD.Stability.CLq =   6.248;
    LD.Stability.Cmq = -18.410;
%deltae --> elevator
    LD.Stability.CDdeltae =  0.005; 
    LD.Stability.CLdeltae =  1.391;
    LD.Stability.Cmdeltae = -2.930;
%deltafr --> rigth flaperon
    LD.Stability.CDdeltafr = LD.Stability.CDdeltae/2;   %TO BE CALCULATED
    LD.Stability.CLdeltafr = LD.Stability.CLdeltae/2;   %TO BE CALCULATED
    LD.Stability.Cmdeltafr = 0.000;           %TO BE CALCULATED
%deltafl --> left flaperon
    LD.Stability.CDdeltafl = LD.Stability.CDdeltae/2;   %TO BE CALCULATED
    LD.Stability.CLdeltafl = LD.Stability.CLdeltae/2;   %TO BE CALCULATED
    LD.Stability.Cmdeltafl = 0.000;           %TO BE CALCULATED
    

%% LATERAL-DIRECTIONAL DERIVATIVES
%beta
    LD.Stability.CYbeta = -0.4250;
    LD.Stability.Clbeta = -0.0286;
    LD.Stability.Cnbeta =  0.0170;
%p
    LD.Stability.CYp = -0.004;
    LD.Stability.Clp = -0.472;
    LD.Stability.Cnp = -0.019;
%r
    LD.Stability.CYr =  0.120;
    LD.Stability.Clr =  0.049;
    LD.Stability.Cnr = -0.040;
%deltar --> rudder
    LD.Stability.CYdeltar = -0.212;
    LD.Stability.Cldeltar = -0.002;
    LD.Stability.Cndeltar =  0.051;
%deltafr --> rigth flaperon
    LD.Stability.CYdeltafr =  0.000;
    LD.Stability.Cldeltafr = -0.125;  %-deltaa/2
    LD.Stability.Cndeltafr =  0.001;  %-deltaa/2
%deltafl --> left flaperon
    LD.Stability.CYdeltafl =  0.000;
    LD.Stability.Cldeltafl =  0.125;  %deltaa/2
    LD.Stability.Cndeltafl = -0.001;  %deltaa/2
        
        
        


