%% CLEAR VARIABLES
clc
close all
clearvars -except LD



%% FLIGHT CONDITION (FC)
FC.hs  =   91.44;                    	  % Altura en [m]
FC.a0  =    340;                          % Velocidad del sonido [m/s]
FC.us  =   42.22;                         % Velocidad condicion de referencia [m/s]
FC.Ms  =   FC.us/FC.a0;      			  % Número de Mach
FC.rho =   1.21;					      % Densidad del aire en [kg/m^3]
FC.qs  =   0.5*FC.rho*(FC.us^2);          % Presion dinámica condicion de referencia [Pa]
FC.gravity   = 9.81;                      % Gravedad en [m/s^2]
 


%% GEOMETRIC DATA (GEO)
Geo.Sw  =  LD.sref;          % Superficie alar [m^2]
Geo.c   =  LD.cbar;          % Cuerda [m]
Geo.b   =  LD.blref;         % envergadura [m]
Geo.alphaT  = deg2rad(0.0); % alpha de los motores (tambien llamado Epsilon)
Geo.alphasb = deg2rad(0.0); % alpha entre ejes estabilidad y cuerpo
Geo.xT = 0; %Distancia de los motores al CdG proyectada en el eje Xs [m] (En este caso no es 0, pero no afecta por ser alphaT=0)
Geo.zT = 0; %Distancia de los motores al CdG proyectada en el eje Zs [m]
Geo.dT = cos(Geo.alphaT)*(Geo.zT/Geo.c) - sin(Geo.alphaT)*(Geo.xT/Geo.c); %Brazo de momentos adimensional de los motores [-]



%% WEIGHT AND BALANCE (WB)
WB.m     = LD.Inertia.mass;       % masa en [kg]
WB.Ixx_b = LD.Inertia.Ixx;        % Ixx en ejes cuerpo [kg*m^2]
WB.Iyy_b = LD.Inertia.Iyy;        % Ixx en ejes cuerpo [kg*m^2]
WB.Izz_b = LD.Inertia.Izz;        % Ixx en ejes cuerpo [kg*m^2]
WB.Jxz_b = LD.Inertia.Ixz;        % Ixx en ejes cuerpo [kg*m^2]

%Transform Inertia moments from body axis to stability axis
Ixxs_Izzs_Jxzs = [     cos(Geo.alphasb)^2,      sin(Geo.alphasb)^2, -sin(2*Geo.alphasb);...
                       sin(Geo.alphasb)^2,      cos(Geo.alphasb)^2,  sin(2*Geo.alphasb);...
                   0.5*sin(2*Geo.alphasb), -0.5*sin(2*Geo.alphasb),  cos(2*Geo.alphasb)] * [WB.Ixx_b; WB.Izz_b; WB.Jxz_b];

WB.Ixx_s = Ixxs_Izzs_Jxzs(1);
WB.Izz_s = Ixxs_Izzs_Jxzs(2);
WB.Jxz_s = Ixxs_Izzs_Jxzs(3);
clear Ixxs_Izzs_Jxzs  



%% REFERENCE CONDITION (RC)
RC.Thetas = 0.000;  %Theta en la condición de referencia
RC.CLs    = WB.m*FC.gravity*cos(RC.Thetas)/(FC.qs*Geo.Sw);   %Coeficiente de sustentacion referencia
RC.CDs    = 0.0276; %Coeficiente de resistencia referencia obtenido de la polar de Raymer
RC.CTs    = 0.0276; %Coeficiente de traccion referencia
RC.Cmas   = 0.000;  %Coeficiente de momento aerodinamico de referencia
RC.CmTs   = 0.000;  %Coeficiente de momentos debido a los motores de referencia



%% ADIMENSIONAL STABILITY DERIVATIVES (ASD) --> LONGITUDINAL
% 0 --> No se usan
    ASD.long.CD_0       =   0.024;
    ASD.long.CL_0       =   0.181;
    ASD.long.Cma_0      =  -0.003;
% u(Solo se tienen en cuenta si volamos en compresible)
    ASD.long.CD_u       =   0.000;
    ASD.long.CL_u       =   0.000;
    ASD.long.Cma_u      =   0.000;
% alpha
    ASD.long.CD_alpha   =   0.239;
    ASD.long.CL_alpha   =   4.138;
    ASD.long.Cma_alpha  =  -0.707;
% alphaDot
    ASD.long.CD_alphaDot  =  0.000;
    ASD.long.CL_alphaDot  =  7.411;  % CL respecto a alpha punto adimensional
    ASD.long.Cma_alphaDot = -4.063;
% q
    ASD.long.CD_q       =   0.000;
    ASD.long.CL_q       =   6.248;
    ASD.long.Cma_q      = -18.410;
% deltae
    ASD.long.CD_deltae  =  0.005;
    ASD.long.CL_deltae  =  1.391;
    ASD.long.Cma_deltae = -2.930;
% deltae_dot   
    ASD.long.CD_deltaeDot  = 0.000;
    ASD.long.CL_deltaeDot  = 0.000;
    ASD.long.Cma_deltaeDot = 0.000;   
% Motores
    ASD.long.CT_u       =   0.000; %To be calculated in detail...
    ASD.long.CT_alpha   =   0.000; %To be calculated in detail...
    ASD.long.CT_deltae  =   0.000; %To be calculated in detail...
    ASD.long.CmT_u      =   0.000;
    ASD.long.CmT_alpha  =   0.000;

    
%Derivadas de Estabilidad ASD longitudinales en ejes estabilidad a partir de las parciales
%(Hipotesis de alpha pequeña)
% Static
    ASD.long.Cxs      =  RC.CTs*cos(Geo.alphasb + Geo.alphaT) - RC.CDs;
    ASD.long.Czs      = -RC.CTs*sin(Geo.alphasb + Geo.alphaT) - RC.CLs;
    ASD.long.Cms      =  RC.Cmas + RC.CmTs;
% u
    ASD.long.Cx_u      =  ASD.long.CT_u * cos(Geo.alphasb + Geo.alphaT) - ASD.long.CD_u;
    ASD.long.Cz_u      = -ASD.long.CL_u - ASD.long.CT_u * sin(Geo.alphasb + Geo.alphaT);
    ASD.long.Cm_u      =  ASD.long.Cma_u + ASD.long.CmT_u;
% alpha
    ASD.long.Cx_alpha  =  ASD.long.CT_alpha*cos(Geo.alphasb + Geo.alphaT) + RC.CLs - ASD.long.CD_alpha;
    ASD.long.Cz_alpha  = -ASD.long.CT_alpha*sin(Geo.alphasb + Geo.alphaT) - ASD.long.CL_alpha - RC.CDs;
    ASD.long.Cm_alpha  =  ASD.long.Cma_alpha + ASD.long.CmT_alpha;
% alphaDot
    ASD.long.Cx_alphaDot = -ASD.long.CD_alphaDot;
    ASD.long.Cz_alphaDot = -ASD.long.CL_alphaDot;
    ASD.long.Cm_alphaDot =  ASD.long.Cma_alphaDot;
% q
    ASD.long.Cx_q      = -ASD.long.CD_q;
    ASD.long.Cz_q      = -ASD.long.CL_q;
    ASD.long.Cm_q      =  ASD.long.Cma_q;
% deltae
    ASD.long.Cx_deltae =  ASD.long.CT_deltae*cos(Geo.alphasb + Geo.alphaT) - ASD.long.CD_deltae;
    ASD.long.Cz_deltae = -ASD.long.CT_deltae*sin(Geo.alphasb + Geo.alphaT) - ASD.long.CL_deltae;
    ASD.long.Cm_deltae =  ASD.long.Cma_deltae;     
% deltae_dot   
    ASD.long.Cx_deltaeDot = -ASD.long.CD_deltaeDot;
    ASD.long.Cz_deltaeDot = -ASD.long.CL_deltaeDot;
    ASD.long.Cm_deltaeDot =  ASD.long.Cma_deltaeDot; 
    


%% ADIMENSIONAL STABILITY DERIVATIVES (ASD) --> LATERAL-DIRECTIONAL
%beta
    ASD.lat.CY_beta   =  -0.425;
    ASD.lat.Cl_beta   =  -0.029;
    ASD.lat.Cn_beta   =   0.017;
%p
    ASD.lat.CY_p      =  -0.004;
    ASD.lat.Cl_p      =  -0.472;
    ASD.lat.Cn_p      =  -0.019;
%r
    ASD.lat.CY_r      =   0.120;
    ASD.lat.Cl_r      =   0.049;
    ASD.lat.Cn_r      =  -0.040;
%deltar
    ASD.lat.CY_deltar =  -0.212;
    ASD.lat.Cl_deltar =  -0.002;
    ASD.lat.Cn_deltar =   0.051;
%deltaa
    ASD.lat.CY_deltaa =   0.000;
    ASD.lat.Cl_deltaa =   0.250;
    ASD.lat.Cn_deltaa =  -0.002;
%deltarDot
    ASD.lat.CY_deltarDot = 0.000;
    ASD.lat.Cl_deltarDot = 0.000;
    ASD.lat.Cn_deltarDot = 0.000;
%deltaaDot
    ASD.lat.CY_deltaaDot = 0.000;
    ASD.lat.Cl_deltaaDot = 0.000;
    ASD.lat.Cn_deltaaDot = 0.000;   
%Motores
    ASD.lat.CT_beta      = 0.000;

    
%% ADIMENSIONAL STABILITY DERIVATIVES (ASD) --> PARAMETERS
   ASD.long.mu = WB.m / (0.5 * FC.rho * Geo.Sw * Geo.c);
   ASD.lat.mu  = WB.m / (0.5 * FC.rho * Geo.Sw * Geo.b);
   ASD.long.Iy = WB.Iyy_b / (FC.rho * Geo.Sw * (Geo.c/2)^3);
   ASD.lat.Ix  = WB.Ixx_s / (FC.rho * Geo.Sw * (Geo.b/2)^3);
   ASD.lat.Iz  = WB.Izz_s / (FC.rho * Geo.Sw * (Geo.b/2)^3);
   ASD.lat.Jxz = WB.Jxz_s / (FC.rho * Geo.Sw * (Geo.b/2)^3);

  
   
%% CHARACTERISTIC EQUATION (ChEq) --> LONGITUDINAL
% (Cuartica de estabilidad longitudinal)
% Metodo adimensional Longitudinal
% Construimos la matriz con las ecuaciones adimensionales hecha su transformada de Laplace
% La primera fila corresponde a "u", la segunda a "alpha" y la tercera a "theta"
% Matriz --> M = MSS*s^2 + MS*s + ML
    syms s

% MSS
    ChEq.long.MSS = zeros (3);
    ChEq.long.MSS(3,3) =  ASD.long.Iy;
    
% MS
    ChEq.long.MS  = zeros (3);
    ChEq.long.MS(1,1) =   2*ASD.long.mu;
    ChEq.long.MS(2,2) =  (2*ASD.long.mu - ASD.long.Cz_alphaDot);
    ChEq.long.MS(2,3) =  -2*ASD.long.mu + ASD.long.Cz_q;
    ChEq.long.MS(3,2) =  -ASD.long.Cm_alphaDot;
    ChEq.long.MS(3,3) =  -ASD.long.Cm_q;

% ML
    ChEq.long.ML  = zeros (3);
    ChEq.long.ML(1,1) =   2*ASD.long.Czs*tan(RC.Thetas) - ASD.long.Cx_u;
    ChEq.long.ML(1,2) = - ASD.long.Cx_alpha;
    ChEq.long.ML(1,3) = - ASD.long.Czs;
    ChEq.long.ML(2,1) = -(2*ASD.long.Czs + ASD.long.Cz_u);
    ChEq.long.ML(2,2) = - ASD.long.Cz_alpha;
    ChEq.long.ML(2,3) = - ASD.long.Czs*tan(RC.Thetas);
    ChEq.long.ML(3,1) = - ASD.long.Cm_u;
    ChEq.long.ML(3,2) = - ASD.long.Cm_alpha;
    ChEq.long.ML(3,3) =   0;

% Matriz --> M = MSS*s^2 + MS*s + ML
% Se obtiene la matriz dimensional en el dominio del tiempo sustituyendo
% por D por t*·d/dt, donde t*=c/(2us) y se dimensionaliza la velocidad
% u = us·u_adim
    ChEq.long.M = ChEq.long.MSS.*(Geo.c/(2*FC.us))^2.*s^2 + ChEq.long.MS.*(Geo.c/(2*FC.us)).*s + ChEq.long.ML;
    
% Calculamos el determinante y los coeficientes de este
% Segun el criterio de Routh, si la cuartica de estabilidad es:
% A*s^4 + B*s^3 + C*s^2 + D*s + E = 0
% El sistema es estable cuando: A>0, B>0, C>0, D>0, E>0 y [D(BC-AD)-(B^2)E]>0
    ChEq.long.det_M  = det(ChEq.long.M);
    ChEq.long.Coeffs = sym2poly(ChEq.long.det_M);
    
% Calculamos los autovalores e identificamos el modo fugoide y corto periodo
    ChEq.long.Eigenvalues = roots(ChEq.long.Coeffs);
    if (length(ChEq.long.Eigenvalues) == 4) 
        for j=1:4; R(j)=isreal(ChEq.long.Eigenvalues(j)); end; clear j %#ok<SAGROW>
        if (sum(R) == 0)
        	ChEq.long.EigenStruct = getEigenData(ChEq.long.Eigenvalues);
            [~,I] = sort(ChEq.long.EigenStruct.freqNat);
            %Phugoid --> Modo Fugoide
            Modes.Phugoid.mode        = ChEq.long.EigenStruct.mode{I(1)};
            Modes.Phugoid.Eigenvalue  =[ChEq.long.EigenStruct.value(I(1)),ChEq.long.EigenStruct.value(I(2))];
            Modes.Phugoid.freqNat     = ChEq.long.EigenStruct.freqNat(I(1));
            Modes.Phugoid.Damp        = ChEq.long.EigenStruct.Damp(I(1));
            Modes.Phugoid.Period      = ChEq.long.EigenStruct.Period(I(1));
            Modes.Phugoid.t_12        = ChEq.long.EigenStruct.t_12(I(1));
            Modes.Phugoid.T2          = ChEq.long.EigenStruct.T2(I(1));
            Modes.Phugoid.Tau         = ChEq.long.EigenStruct.Tau(I(1));
            syms u alpha theta
            eq = double(subs(ChEq.long.M,s,ChEq.long.EigenStruct.value(I(2))))*[u alpha theta]';
            eq = subs(eq,theta,1);
            sol = solve([eq(1),eq(2)],[u,alpha]);
            Modes.Phugoid.Eigenvector = [double(sol.u) double(sol.alpha) 1];
            

            %Short Period --> Corto Periodo
            Modes.ShortPeriod.mode  = ChEq.long.EigenStruct.mode{I(3)};
            Modes.ShortPeriod.value = [ChEq.long.EigenStruct.value(I(3)),ChEq.long.EigenStruct.value(I(4))];
            Modes.ShortPeriod.freqNat = ChEq.long.EigenStruct.freqNat(I(3));
            Modes.ShortPeriod.Damp = ChEq.long.EigenStruct.Damp(I(3));
            Modes.ShortPeriod.Period = ChEq.long.EigenStruct.Period(I(3));
            Modes.ShortPeriod.t_12 = ChEq.long.EigenStruct.t_12(I(3));
            Modes.ShortPeriod.T2 = ChEq.long.EigenStruct.T2(I(3));
            Modes.ShortPeriod.Tau = ChEq.long.EigenStruct.Tau(I(3));
            eq = double(subs(ChEq.long.M,s,ChEq.long.EigenStruct.value(I(4))))*[u alpha theta]';
            eq = subs(eq,theta,1);
            sol = solve([eq(1),eq(2)],[u,alpha]);
            Modes.ShortPeriod.Eigenvector = [double(sol.u) double(sol.alpha) 1];
            clear I u alpha theta eq sol
        else
            wrn = msgbox({'Alguno de los modos longitudinales no es oscilatorio.',...
                          'Se recomienda revisar a mano lo que está sucediendo.'},...
                          'Warning','warn');
            uiwait(wrn);
            disp('Alguno de los modos longitudinales no es oscilatorio.')
            disp('Se recomienda revisar a mano lo que está sucediendo.')
            pause
        end
        clear R
    else
        wrn = msgbox({'El número de modos longitudinales es diferente de 4.',...
                      'Se recomienda revisar a mano lo que está sucediendo.'},...
                      'Warning','warn');
        uiwait(wrn);
        disp('El número de modos longitudinales es diferente de 4.')
        disp('Se recomienda revisar a mano lo que está sucediendo.')
        pause
    end
 
    
    
    
%% CHARACTERISTIC EQUATION (ChEq) --> LATERAL-DIRECTIONAL
% (Cuartica de estabilidad lateral-direccional)
% Metodo adimensional Lateral-Direccional
% Construimos la matriz con las ecuaciones adimensionales hecha su transformada de Laplace
% La primera fila corresponde a "beta", la segunda a "Phi" y la tercera a "Psi"
% Matriz --> M = MSS*s^2 + MS*s + ML
    syms s

% MSS
    ChEq.lat.MSS = zeros (3);
    ChEq.lat.MSS(2,2) =  ASD.lat.Ix;
    ChEq.lat.MSS(2,3) = -ASD.lat.Jxz;
    ChEq.lat.MSS(3,2) = -ASD.lat.Jxz;
    ChEq.lat.MSS(3,3) =  ASD.lat.Iz;
    
% MS
    ChEq.lat.MS  = zeros (3);
    ChEq.lat.MS(1,1) =   2*ASD.lat.mu;
    ChEq.lat.MS(1,2) = - ASD.lat.CY_p;
    ChEq.lat.MS(1,3) =   2*ASD.lat.mu - ASD.lat.CY_r;
    ChEq.lat.MS(2,2) = - ASD.lat.Cl_p;
    ChEq.lat.MS(2,3) = - ASD.lat.Cl_r;
    ChEq.lat.MS(3,2) = - ASD.lat.Cn_p;
    ChEq.lat.MS(3,3) = - ASD.lat.Cn_r;

% ML
    ChEq.lat.ML  = zeros (3);
    ChEq.lat.ML(1,1) = - ASD.lat.CY_beta;
    ChEq.lat.ML(1,2) = - RC.CLs;
    ChEq.lat.ML(2,1) = - ASD.lat.Cl_beta;
    ChEq.lat.ML(3,1) = - ASD.lat.Cn_beta;

    
% Matriz --> M = MSS*s^2 + MS*s + ML
% Se obtiene la matriz dimensional en el dominio del tiempo sustituyendo
% por D por t*·d/dt, donde t*=b/(2us) y se dimensionaliza la velocidad
% u = us·u_adim
    ChEq.lat.M = ChEq.lat.MSS.*(Geo.b/(2*FC.us))^2.*s^2 + ChEq.lat.MS.*(Geo.b/(2*FC.us)).*s + ChEq.lat.ML;

% Calculamos el determinante y los coeficientes de este
% Segun el criterio de Routh, si la cuartica de estabilidad es:
% A*s^4 + B*s^3 + C*s^2 + D*s + E = 0
% El sistema es estable cuando: A>0, B>0, C>0, D>0, E>0 y [D(BC-AD)-(B^2)E]>0
    ChEq.lat.det_M  = det(ChEq.lat.M);
    ChEq.lat.Coeffs = sym2poly(ChEq.lat.det_M); 
    
% Calculamos los autovalores e identificamos los modos espira, convergencia
% en balance y balanceo holandes
    ChEq.lat.Eigenvalues = roots(ChEq.lat.Coeffs);
    
    if (length(ChEq.lat.Eigenvalues) == 5) 
        for j=1:5; R(j)=isreal(ChEq.lat.Eigenvalues(j)); end; clear j %#ok<SAGROW>
        if (sum(R) == 3)
        	ChEq.lat.EigenStruct = getEigenData(ChEq.lat.Eigenvalues);
            [~,I] = sort(abs(real(ChEq.lat.EigenStruct.value)));
            %Spiral --> Modo Espiral
            Modes.Spiral.mode  = ChEq.lat.EigenStruct.mode{I(2)};
            Modes.Spiral.value = ChEq.lat.EigenStruct.value(I(2));
            Modes.Spiral.freqNat = ChEq.lat.EigenStruct.freqNat(I(2));
            Modes.Spiral.Damp = ChEq.lat.EigenStruct.Damp(I(2));
            Modes.Spiral.Period = ChEq.lat.EigenStruct.Period(I(2));
            Modes.Spiral.t_12 = ChEq.lat.EigenStruct.t_12(I(2));
            Modes.Spiral.T2 = ChEq.lat.EigenStruct.T2(I(2));
            Modes.Spiral.Tau = ChEq.lat.EigenStruct.Tau(I(2));
            syms beta Phi Psi
            eq = double(subs(ChEq.lat.M,s,ChEq.lat.EigenStruct.value(I(2))))*[beta Phi Psi]';
            eq = subs(eq,Phi,1);
            sol = solve([eq(1),eq(3)],[beta Psi]);
            Modes.Spiral.Eigenvector = [double(sol.beta) 1 double(sol.Psi)];
            
            %Roll --> Modo de Convergencia en Balance
            Modes.RollingConvergence.mode  = ChEq.lat.EigenStruct.mode{I(5)};
            Modes.RollingConvergence.value = ChEq.lat.EigenStruct.value(I(5));
            Modes.RollingConvergence.freqNat = ChEq.lat.EigenStruct.freqNat(I(5));
            Modes.RollingConvergence.Damp = ChEq.lat.EigenStruct.Damp(I(5));
            Modes.RollingConvergence.Period = ChEq.lat.EigenStruct.Period(I(5));
            Modes.RollingConvergence.t_12 = ChEq.lat.EigenStruct.t_12(I(5));
            Modes.RollingConvergence.T2 = ChEq.lat.EigenStruct.T2(I(5));
            Modes.RollingConvergence.Tau = ChEq.lat.EigenStruct.Tau(I(5));
            eq = double(subs(ChEq.lat.M,s,ChEq.lat.EigenStruct.value(I(5))))*[beta Phi Psi]';
            eq = subs(eq,Phi,1);
            sol = solve([eq(1),eq(3)],[beta Psi]);
            Modes.RollingConvergence.Eigenvector = [double(sol.beta) 1 double(sol.Psi)];
            
            %Dutch Roll --> Modo de Balanceo Holandes
            Modes.DutchRoll.mode  = ChEq.lat.EigenStruct.mode{I(3)};
            Modes.DutchRoll.value = [ChEq.lat.EigenStruct.value(I(3)),ChEq.long.EigenStruct.value(I(4))];
            Modes.DutchRoll.freqNat = ChEq.lat.EigenStruct.freqNat(I(3));
            Modes.DutchRoll.Damp = ChEq.lat.EigenStruct.Damp(I(3));
            Modes.DutchRoll.Period = ChEq.lat.EigenStruct.Period(I(3));
            Modes.DutchRoll.t_12 = ChEq.lat.EigenStruct.t_12(I(3));
            Modes.DutchRoll.T2 = ChEq.lat.EigenStruct.T2(I(3));
            Modes.DutchRoll.Tau = ChEq.lat.EigenStruct.Tau(I(3));
            eq = double(subs(ChEq.lat.M,s,ChEq.lat.EigenStruct.value(I(4))))*[beta Phi Psi]';
            eq = subs(eq,Phi,1);
            sol = solve([eq(1),eq(3)],[beta Psi]);
            Modes.DutchRoll.Eigenvector = [double(sol.beta) 1 double(sol.Psi)];
            clear I beta Phi Psi eq sol
        else
            wrn = msgbox({'Los modos lateral-direccionales no son 3 reales y 2 oscilatorios.',...
                          'Se recomienda revisar a mano lo que está sucediendo.'},...
                          'Warning','warn');
            uiwait(wrn);
            disp('Los modos lateral-direccionales no son 3 reales y 2 oscilatorios.')
            disp('Se recomienda revisar a mano lo que está sucediendo.')
            pause
        end
        clear R
    else
        wrn = msgbox({'El número de modos lateral-direccionales es diferente de 5.',...
                      'Se recomienda revisar a mano lo que está sucediendo.'},...
                      'Warning','warn');
        uiwait(wrn);
        disp('El número de modos lateral-direccionales es diferente de 5.')
        disp('Se recomienda revisar a mano lo que está sucediendo.')
        pause
    end
    
    
    
%% TRANSFER FUNCTIONS (TF) --> LONGITUDINAL

% FUNCIÓN DE TRANSFERENCIA DE ELEVADOR A VELOCIDAD
    % Construimos el numerador
    % MK2
        TF.long.MK2 = ChEq.long.MSS;
        TF.long.MK2(:,1) = 0 ;
    % MK1
        TF.long.MK1 = ChEq.long.MS;
        TF.long.MK1(:,1)  = 0 ;
        TF.long.MK1(1,1)  = ASD.long.Cx_deltaeDot;
        TF.long.MK1(2,1)  = ASD.long.Cz_deltaeDot;
        TF.long.MK1(3,1)  = ASD.long.Cm_deltaeDot;
    % MK0
        TF.long.MK0 = ChEq.long.ML;
        TF.long.MK0(:,1)   = 0 ;
        TF.long.MK0(1,1)   = ASD.long.Cx_deltae;
        TF.long.MK0(2,1)   = ASD.long.Cz_deltae;
        TF.long.MK0(3,1)   = ASD.long.Cm_deltae;
    % Matriz --> MK = MK2*s^2 + MK1*s + MK0
    % Se obtiene la matriz dimensional en el dominio del tiempo sustituyendo
    % por D por t*·d/dt, donde t*=c/(2us)
        TF.long.MK = TF.long.MK2.*(Geo.c/(2*FC.us)).^2.*s^2 + TF.long.MK1.*(Geo.c/(2*FC.us)).*s + TF.long.MK0;
    % Calculamos el determinante y los coeficientes de este
        TF.long.det_MK = det(TF.long.MK);
        TF.long.Coeffs = sym2poly(TF.long.det_MK); 
    % Se obtiene la función de transferencia de elevador a velocidad adimensional 
        TF.long.Elevator_Speed = tf(TF.long.Coeffs,ChEq.long.Coeffs);
    % Se dimensionaliza la funcion de transferencia para ser elevador a velocidad
        TF.long.Elevator_Speed = TF.long.Elevator_Speed * FC.us;
        disp('Speed to Elevator transfer function [(m/s)/rad]:')
        TF.long.Elevator_Speed
        fprintf('\n\n')

        
        
        
% FUNCIÓN DE TRANSFERENCIA DE ELEVADOR A ANGULO DE ATAQUE
    % Construimos el numerador
    % MK2
        TF.long.MK2 = ChEq.long.MSS;
        TF.long.MK2(:,2) = 0 ;
    % MK1
        TF.long.MK1 = ChEq.long.MS;
        TF.long.MK1(:,2)  = 0 ;
        TF.long.MK1(1,2)  = ASD.long.Cx_deltaeDot;
        TF.long.MK1(2,2)  = ASD.long.Cz_deltaeDot;
        TF.long.MK1(3,2)  = ASD.long.Cm_deltaeDot;
    % MK0
        TF.long.MK0 = ChEq.long.ML;
        TF.long.MK0(:,2)   = 0 ;
        TF.long.MK0(1,2)   = ASD.long.Cx_deltae;
        TF.long.MK0(2,2)   = ASD.long.Cz_deltae;
        TF.long.MK0(3,2)   = ASD.long.Cm_deltae;
    % Matriz --> MK = MK2*s^2 + MK1*s + MK0
    % Se obtiene la matriz dimensional en el dominio del tiempo sustituyendo
    % por D por t*·d/dt, donde t*=c/(2us)
        TF.long.MK = TF.long.MK2.*(Geo.c/(2*FC.us)).^2.*s^2 + TF.long.MK1.*(Geo.c/(2*FC.us)).*s + TF.long.MK0;
    % Calculamos el determinante y los coeficientes de este
        TF.long.det_MK = det(TF.long.MK);
        TF.long.Coeffs = sym2poly(TF.long.det_MK); 
    % Se obtiene la función de transferencia de elevador a angulo de ataque    
        TF.long.Elevator_Alpha = tf(TF.long.Coeffs,ChEq.long.Coeffs);
        disp('Alpha to Elevator transfer function [rad/rad]:')
        TF.long.Elevator_Alpha
        fprintf('\n\n')

        
        
        
 % FUNCIÓN DE TRANSFERENCIA DE ELEVADOR A ANGULO DE ASIENTO
    % Construimos el numerador
    % MK2
        TF.long.MK2 = ChEq.long.MSS;
        TF.long.MK2(:,3) = 0 ;
    % MK1
        TF.long.MK1 = ChEq.long.MS;
        TF.long.MK1(:,3)  = 0 ;
        TF.long.MK1(1,3)  = ASD.long.Cx_deltaeDot;
        TF.long.MK1(2,3)  = ASD.long.Cz_deltaeDot;
        TF.long.MK1(3,3)  = ASD.long.Cm_deltaeDot;
    % MK0
        TF.long.MK0 = ChEq.long.ML;
        TF.long.MK0(:,3)   = 0 ;
        TF.long.MK0(1,3)   = ASD.long.Cx_deltae;
        TF.long.MK0(2,3)   = ASD.long.Cz_deltae;
        TF.long.MK0(3,3)   = ASD.long.Cm_deltae;
    % Matriz --> MK = MK2*s^2 + MK1*s + MK0
    % Se obtiene la matriz dimensional en el dominio del tiempo sustituyendo
    % por D por t*·d/dt, donde t*=c/(2us)
        TF.long.MK = TF.long.MK2.*(Geo.c/(2*FC.us)).^2.*s^2 + TF.long.MK1.*(Geo.c/(2*FC.us)).*s + TF.long.MK0;
    % Calculamos el determinante y los coeficientes de este
        TF.long.det_MK = det(TF.long.MK);
        TF.long.Coeffs = sym2poly(TF.long.det_MK); 
    % Se obtiene la función de transferencia de elevador a angulo de asiento    
        TF.long.Elevator_Theta = tf(TF.long.Coeffs,ChEq.long.Coeffs);
        disp('Theta to Elevator transfer function [rad/rad]:')
        TF.long.Elevator_Theta
        fprintf('\n\n')
        
        
        
        
 % FUNCIÓN DE TRANSFERENCIA DE ELEVADOR A VELOCIDAD ANGULAR DE CABECEO (q)
    % La funcion de transferencia de elevador a velocidad angular de
    % cabeceo adimensional es la función de transferecia de elevador a
    % ángulo de asiento multiplicada por la variable de Laplace (s)
        TF.long.Elevator_ThetaRate = tf([TF.long.Coeffs,0],ChEq.long.Coeffs);
    % Devolvemos las dimensiones a la velocidad angular de cabeceo, para
    % ello: q = q_adim * (2·Us/c)
        TF.long.Elevator_ThetaRate = TF.long.Elevator_ThetaRate * (2*FC.us/Geo.c);
        disp('Theta Rate to Elevator transfer function [(rad/s)/rad]:')
        TF.long.Elevator_ThetaRate
        fprintf('\n\n')
        
        
        
    
% Clear unnecessary fields
	TF.long = rmfield(TF.long,{'MK', 'MK0', 'MK1', 'MK2', 'det_MK', 'Coeffs'});
    
    
        
% SE CALCULAN LAS GANANCIAS ESTATICAS
%Elevador
    response = step(TF.long.Elevator_Speed);
    TF.long.StaticGains.Elevator_Speed = response(end);
    
    response = step(TF.long.Elevator_Alpha);
    TF.long.StaticGains.Elevator_Alpha = response(end);
    
    response = step(TF.long.Elevator_Theta);
    TF.long.StaticGains.Elevator_Theta = response(end);
    
    response = step(TF.long.Elevator_ThetaRate);
    TF.long.StaticGains.Elevator_ThetaRate = response(end);
    clear response
    
    
% SE GUARDAN LOS NUMERADORES Y EL DENOMINADOR (Por comodidad)
    TF.long.Numerators  = [TF.long.Elevator_Speed.Numerator,...
                           TF.long.Elevator_Alpha.Numerator,...
                           TF.long.Elevator_Theta.Numerator,...
                           TF.long.Elevator_ThetaRate.Numerator];
    TF.long.Denominator = ChEq.long.Coeffs;

    

    
    
%% TRANSFER FUNCTIONS (TF) --> LATERAL-DIRECTIONAL

% FUNCIÓN DE TRANSFERENCIA DE ALERONES A ÁNGULO DE RESBALAMIENTO
    % Construimos el numerador
    % MK2
        TF.lat.MK2 = ChEq.lat.MSS;
        TF.lat.MK2(:,1) = 0 ;
    % MK1
        TF.lat.MK1 = ChEq.lat.MS;
        TF.lat.MK1(:,1)  = 0 ;
        TF.lat.MK1(1,1)  = ASD.lat.CY_deltaaDot;
        TF.lat.MK1(2,1)  = ASD.lat.Cl_deltaaDot;
        TF.lat.MK1(3,1)  = ASD.lat.Cn_deltaaDot;
    % MK0
        TF.lat.MK0 = ChEq.lat.ML;
        TF.lat.MK0(:,1)   = 0 ;
        TF.lat.MK0(1,1)   = ASD.lat.CY_deltaa;
        TF.lat.MK0(2,1)   = ASD.lat.Cl_deltaa;
        TF.lat.MK0(3,1)   = ASD.lat.Cn_deltaa;
    % Matriz --> MK = MK2*s^2 + MK1*s + MK0
    % Se obtiene la matriz dimensional en el dominio del tiempo sustituyendo
    % por D por t*·d/dt, donde t*=c/(2us)
        TF.lat.MK = TF.lat.MK2.*(Geo.b/(2*FC.us)).^2.*s^2 + TF.lat.MK1.*(Geo.b/(2*FC.us)).*s + TF.lat.MK0;
    % Calculamos el determinante y los coeficientes de este
        TF.lat.det_MK = det(TF.lat.MK);
        TF.lat.Coeffs = sym2poly(TF.lat.det_MK); 
    % Se obtiene la función de transferencia de alerones a angulo de resbalamiento 
        TF.lat.Aileron_Sideslip = tf(TF.lat.Coeffs,ChEq.lat.Coeffs);
        disp('Sideslip to Aileron transfer function [rad/rad]:')
        TF.lat.Aileron_Sideslip
        fprintf('\n\n')

        
        
        
% FUNCIÓN DE TRANSFERENCIA DE ALERONES A ÁNGULO DE BALANCE
    % Construimos el numerador
    % MK2
        TF.lat.MK2 = ChEq.lat.MSS;
        TF.lat.MK2(:,2) = 0 ;
    % MK1
        TF.lat.MK1 = ChEq.lat.MS;
        TF.lat.MK1(:,2)  = 0 ;
        TF.lat.MK1(1,2)  = ASD.lat.CY_deltaaDot;
        TF.lat.MK1(2,2)  = ASD.lat.Cl_deltaaDot;
        TF.lat.MK1(3,2)  = ASD.lat.Cn_deltaaDot;
    % MK0
        TF.lat.MK0 = ChEq.lat.ML;
        TF.lat.MK0(:,2)   = 0 ;
        TF.lat.MK0(1,2)   = ASD.lat.CY_deltaa;
        TF.lat.MK0(2,2)   = ASD.lat.Cl_deltaa;
        TF.lat.MK0(3,2)   = ASD.lat.Cn_deltaa;
    % Matriz --> MK = MK2*s^2 + MK1*s + MK0
    % Se obtiene la matriz dimensional en el dominio del tiempo sustituyendo
    % por D por t*·d/dt, donde t*=c/(2us)
        TF.lat.MK = TF.lat.MK2.*(Geo.b/(2*FC.us)).^2.*s^2 + TF.lat.MK1.*(Geo.b/(2*FC.us)).*s + TF.lat.MK0;
    % Calculamos el determinante y los coeficientes de este
        TF.lat.det_MK = det(TF.lat.MK);
        TF.lat.Coeffs = sym2poly(TF.lat.det_MK); 
    % Se obtiene la función de transferencia de alerones a angulo de balance 
        TF.lat.Aileron_Phi = tf(TF.lat.Coeffs,ChEq.lat.Coeffs);
        disp('Phi to Aileron transfer function [rad/rad]:')
        TF.lat.Aileron_Phi
        fprintf('\n\n')

        
        
        
% FUNCIÓN DE TRANSFERENCIA DE ALERONES A VELOCIDAD ANGULAR DE BALANCE
    % La funcion de transferencia de alerones a velocidad angular de
    % balance adimensional es la función de transferecia de alerones a
    % ángulo de balance multiplicada por la variable de Laplace (s)
        TF.lat.Aileron_RollRate = tf([TF.lat.Coeffs,0],ChEq.lat.Coeffs);
    % Devolvemos las dimensiones a la velocidad angular de balance, para
    % ello: p = p_adim * (2·Us/b)
        TF.lat.Aileron_RollRate = TF.lat.Aileron_RollRate * (2*FC.us/Geo.b);
        disp('Roll Rate to Aileron transfer function [(rad/s)/rad]:')
        TF.lat.Aileron_RollRate
        fprintf('\n\n')
        

        
        
% FUNCIÓN DE TRANSFERENCIA DE ALERONES A GUIÑADA
    % Construimos el numerador
    % MK2
        TF.lat.MK2 = ChEq.lat.MSS;
        TF.lat.MK2(:,3) = 0 ;
    % MK1
        TF.lat.MK1 = ChEq.lat.MS;
        TF.lat.MK1(:,3)  = 0 ;
        TF.lat.MK1(1,3)  = ASD.lat.CY_deltaaDot;
        TF.lat.MK1(2,3)  = ASD.lat.Cl_deltaaDot;
        TF.lat.MK1(3,3)  = ASD.lat.Cn_deltaaDot;
    % MK0
        TF.lat.MK0 = ChEq.lat.ML;
        TF.lat.MK0(:,3)   = 0 ;
        TF.lat.MK0(1,3)   = ASD.lat.CY_deltaa;
        TF.lat.MK0(2,3)   = ASD.lat.Cl_deltaa;
        TF.lat.MK0(3,3)   = ASD.lat.Cn_deltaa;
    % Matriz --> MK = MK2*s^2 + MK1*s + MK0
    % Se obtiene la matriz dimensional en el dominio del tiempo sustituyendo
    % por D por t*·d/dt, donde t*=c/(2us)
        TF.lat.MK = TF.lat.MK2.*(Geo.b/(2*FC.us)).^2.*s^2 + TF.lat.MK1.*(Geo.b/(2*FC.us)).*s + TF.lat.MK0;
    % Calculamos el determinante y los coeficientes de este
        TF.lat.det_MK = det(TF.lat.MK);
        TF.lat.Coeffs = sym2poly(TF.lat.det_MK); 
    % Se obtiene la función de transferencia de alerones a angulo de guiñada
        TF.lat.Aileron_Yaw = tf(TF.lat.Coeffs,ChEq.lat.Coeffs);
        disp('Yaw to Aileron transfer function [rad/rad]:')
        TF.lat.Aileron_Yaw
        fprintf('\n\n')

        
        
        
% FUNCIÓN DE TRANSFERENCIA DE ALERONES A VELOCIDAD ANGULAR DE GUIÑADA
    % La funcion de transferencia de alerones a ángulo de guiñada es la
    % función de transferecia de alerones a velocidad angular de guiñada
    % adimensional dividida por la variable de Laplace (s)
        TF.lat.Aileron_YawRate = tf([TF.lat.Coeffs,0],ChEq.lat.Coeffs);
    % Se dimensionaliza la funcion de transferencia para ser alerones a velocidad angular
        TF.lat.Aileron_YawRate = TF.lat.Aileron_YawRate * (2*FC.us/Geo.b);
        disp('Yaw Rate to Aileron transfer function [(rad/s)/rad]:')
        TF.lat.Aileron_YawRate
        fprintf('\n\n')

        
        
        
        
% FUNCIÓN DE TRANSFERENCIA DE TIMÓN DE DIRECCIÓN A ÁNGULO DE RESBALAMIENTO
    % Construimos el numerador
    % MK2
        TF.lat.MK2 = ChEq.lat.MSS;
        TF.lat.MK2(:,1) = 0 ;
    % MK1
        TF.lat.MK1 = ChEq.lat.MS;
        TF.lat.MK1(:,1)  = 0 ;
        TF.lat.MK1(1,1)  = ASD.lat.CY_deltarDot;
        TF.lat.MK1(2,1)  = ASD.lat.Cl_deltarDot;
        TF.lat.MK1(3,1)  = ASD.lat.Cn_deltarDot;
    % MK0
        TF.lat.MK0 = ChEq.lat.ML;
        TF.lat.MK0(:,1)   = 0 ;
        TF.lat.MK0(1,1)   = ASD.lat.CY_deltar;
        TF.lat.MK0(2,1)   = ASD.lat.Cl_deltar;
        TF.lat.MK0(3,1)   = ASD.lat.Cn_deltar;
    % Matriz --> MK = MK2*s^2 + MK1*s + MK0
    % Se obtiene la matriz dimensional en el dominio del tiempo sustituyendo
    % por D por t*·d/dt, donde t*=c/(2us)
        TF.lat.MK = TF.lat.MK2.*(Geo.b/(2*FC.us)).^2.*s^2 + TF.lat.MK1.*(Geo.b/(2*FC.us)).*s + TF.lat.MK0;
    % Calculamos el determinante y los coeficientes de este
        TF.lat.det_MK = det(TF.lat.MK);
        TF.lat.Coeffs = sym2poly(TF.lat.det_MK); 
    % Se obtiene la función de transferencia de timón de dirección a angulo de resbalamiento 
        TF.lat.Rudder_Sideslip = tf(TF.lat.Coeffs,ChEq.lat.Coeffs);
        disp('Sideslip to Rudder transfer function [rad/rad]:')
        TF.lat.Rudder_Sideslip
        fprintf('\n\n')

        
        
        
% FUNCIÓN DE TRANSFERENCIA DE TIMÓN DE DIRECCIÓN A ÁNGULO DE BALANCE
    % Construimos el numerador
    % MK2
        TF.lat.MK2 = ChEq.lat.MSS;
        TF.lat.MK2(:,2) = 0 ;
    % MK1
        TF.lat.MK1 = ChEq.lat.MS;
        TF.lat.MK1(:,2)  = 0 ;
        TF.lat.MK1(1,2)  = ASD.lat.CY_deltarDot;
        TF.lat.MK1(2,2)  = ASD.lat.Cl_deltarDot;
        TF.lat.MK1(3,2)  = ASD.lat.Cn_deltarDot;
    % MK0
        TF.lat.MK0 = ChEq.lat.ML;
        TF.lat.MK0(:,2)   = 0 ;
        TF.lat.MK0(1,2)   = ASD.lat.CY_deltar;
        TF.lat.MK0(2,2)   = ASD.lat.Cl_deltar;
        TF.lat.MK0(3,2)   = ASD.lat.Cn_deltar;
    % Matriz --> MK = MK2*s^2 + MK1*s + MK0
    % Se obtiene la matriz dimensional en el dominio del tiempo sustituyendo
    % por D por t*·d/dt, donde t*=c/(2us)
        TF.lat.MK = TF.lat.MK2.*(Geo.b/(2*FC.us)).^2.*s^2 + TF.lat.MK1.*(Geo.b/(2*FC.us)).*s + TF.lat.MK0;
    % Calculamos el determinante y los coeficientes de este
        TF.lat.det_MK = det(TF.lat.MK);
        TF.lat.Coeffs = sym2poly(TF.lat.det_MK); 
    % Se obtiene la función de transferencia de timón de dirección a angulo de balance 
        TF.lat.Rudder_Phi = tf(TF.lat.Coeffs,ChEq.lat.Coeffs);
        disp('Phi to Rudder transfer function [rad/rad]:')
        TF.lat.Rudder_Phi
        fprintf('\n\n')

        
        
        
% FUNCIÓN DE TRANSFERENCIA DE TIMÓN DE DIRECCIÓN A VELOCIDAD ANGULAR DE BALANCE
    % La funcion de transferencia de timón de dirección a velocidad angular de
    % balance adimensional es la función de transferecia de timón de dirección a
    % ángulo de balance multiplicada por la variable de Laplace (s)
        TF.lat.Rudder_RollRate = tf([TF.lat.Coeffs,0],ChEq.lat.Coeffs);
    % Devolvemos las dimensiones a la velocidad angular de balance, para
    % ello: p = p_adim * (2·Us/b)
        TF.lat.Rudder_RollRate = TF.lat.Rudder_RollRate * (2*FC.us/Geo.b);
        disp('Roll Rate to Rudder transfer function [(rad/s)/rad]:')
        TF.lat.Rudder_RollRate
        fprintf('\n\n')
        

        
        
% FUNCIÓN DE TRANSFERENCIA DE TIMÓN DE DIRECCIÓN A GUIÑADA
    % Construimos el numerador
    % MK2
        TF.lat.MK2 = ChEq.lat.MSS;
        TF.lat.MK2(:,3) = 0 ;
    % MK1
        TF.lat.MK1 = ChEq.lat.MS;
        TF.lat.MK1(:,3)  = 0 ;
        TF.lat.MK1(1,3)  = ASD.lat.CY_deltarDot;
        TF.lat.MK1(2,3)  = ASD.lat.Cl_deltarDot;
        TF.lat.MK1(3,3)  = ASD.lat.Cn_deltarDot;
    % MK0
        TF.lat.MK0 = ChEq.lat.ML;
        TF.lat.MK0(:,3)   = 0 ;
        TF.lat.MK0(1,3)   = ASD.lat.CY_deltar;
        TF.lat.MK0(2,3)   = ASD.lat.Cl_deltar;
        TF.lat.MK0(3,3)   = ASD.lat.Cn_deltar;
    % Matriz --> MK = MK2*s^2 + MK1*s + MK0
    % Se obtiene la matriz dimensional en el dominio del tiempo sustituyendo
    % por D por t*·d/dt, donde t*=c/(2us)
        TF.lat.MK = TF.lat.MK2.*(Geo.b/(2*FC.us)).^2.*s^2 + TF.lat.MK1.*(Geo.b/(2*FC.us)).*s + TF.lat.MK0;
    % Calculamos el determinante y los coeficientes de este
        TF.lat.det_MK = det(TF.lat.MK);
        TF.lat.Coeffs = sym2poly(TF.lat.det_MK); 
    % Se obtiene la función de transferencia de timón de dirección a angulo de guiñada
        TF.lat.Rudder_Yaw = tf(TF.lat.Coeffs,ChEq.lat.Coeffs);
        disp('Yaw to Rudder transfer function [rad/rad]:')
        TF.lat.Rudder_Yaw
        fprintf('\n\n')

        
        
        
% FUNCIÓN DE TRANSFERENCIA DE TIMÓN DE DIRECCIÓN A VELOCIDAD ANGULAR DE GUIÑADA
    % La funcion de transferencia de timón de dirección a ángulo de guiñada es la
    % función de transferecia de timón de dirección a velocidad angular de guiñada
    % adimensional dividida por la variable de Laplace (s)
        TF.lat.Rudder_YawRate = tf([TF.lat.Coeffs,0],ChEq.lat.Coeffs);
    % Se dimensionaliza la funcion de transferencia para ser timón de dirección a velocidad angular
        TF.lat.Rudder_YawRate = TF.lat.Rudder_YawRate * (2*FC.us/Geo.b);
        disp('Yaw Rate to Rudder transfer function [(rad/s)/rad]:')
        TF.lat.Rudder_YawRate
        fprintf('\n\n')        
        
        
        
        
% Clear unnecessary fields
	TF.lat = rmfield(TF.lat,{'MK', 'MK0', 'MK1', 'MK2', 'det_MK', 'Coeffs'});         
        
    
    
% SE CALCULAN LAS GANANCIAS ESTATICAS
% Alerones
    response = step(TF.lat.Aileron_Sideslip);
    TF.lat.StaticGains.Aileron_Sideslip = response(end);
    
    response = step(TF.lat.Aileron_Phi);
    TF.lat.StaticGains.Aileron_Phi = response(end);
    
    response = step(TF.lat.Aileron_RollRate);
    TF.lat.StaticGains.Aileron_RollRate = response(end);
    
    response = step(TF.lat.Aileron_Yaw);
    TF.lat.StaticGains.Aileron_Yaw = response(end);
    
    response = step(TF.lat.Aileron_YawRate);
    TF.lat.StaticGains.Aileron_YawRate = response(end);
    
% Rudder
    response = step(TF.lat.Rudder_Sideslip);
    TF.lat.StaticGains.Rudder_Sideslip = response(end);
    
    response = step(TF.lat.Rudder_Phi);
    TF.lat.StaticGains.Rudder_Phi = response(end);
    
    response = step(TF.lat.Rudder_RollRate);
    TF.lat.StaticGains.Rudder_RollRate = response(end);
    
    response = step(TF.lat.Rudder_Yaw);
    TF.lat.StaticGains.Rudder_Yaw  = response(end);
    
    response = step(TF.lat.Rudder_YawRate);
    TF.lat.StaticGains.Rudder_YawRate = response(end);    
    clear response
    
    
% SE GUARDAN LOS NUMERADORES Y EL DENOMINADOR (Por comodidad)
    TF.lat.Aileron_Numerators  = [TF.lat.Aileron_Sideslip.Numerator,...
                                  TF.lat.Aileron_Phi.Numerator,...
                                  TF.lat.Aileron_RollRate.Numerator,...
                                  TF.lat.Aileron_Yaw.Numerator,...
                                  TF.lat.Aileron_YawRate.Numerator];
                              
    TF.lat.Rudder_Numerators = [TF.lat.Rudder_Sideslip.Numerator,...
                                TF.lat.Rudder_Phi.Numerator,...
                                TF.lat.Rudder_RollRate.Numerator,...
                                TF.lat.Rudder_Yaw.Numerator,...
                                TF.lat.Rudder_YawRate.Numerator];
                            
    TF.lat.Denominator = ChEq.lat.Coeffs;       
    
 % Se guardan las funciones de transferencia en un archivo para ser
 % utilizadas de manera independiente
 
    LibisTF = TF;
    save('.\utilities\LibisTransferFunctions.mat', 'LibisTF')
    
     
%% OTHER USEFUL FUNCTIONS USED IN THE SCRIPT
function [ eigenStruct ] = getEigenData( eigenvalue )
    %GETEIGENDATA Summary of this function goes here
    %   Detailed explanation goes here

    if length(eigenvalue)~=1
        eigenStruct.mode = cell(1,length(eigenvalue));
        eigenStruct.value(1:length(eigenvalue))   = NaN;
        eigenStruct.t_12(1:length(eigenvalue))    = NaN;
        eigenStruct.T2(1:length(eigenvalue))      = NaN;
        eigenStruct.Period(1:length(eigenvalue))  = NaN;
        eigenStruct.freqNat(1:length(eigenvalue)) = NaN;
        eigenStruct.Damp(1:length(eigenvalue))    = NaN;
        eigenStruct.Tau(1:length(eigenvalue))     = NaN;
    end

    for i=1:length(eigenvalue)

        if isreal(eigenvalue(i))
            eigenStruct.mode{i}  = 'Exponential';
            eigenStruct.value(i) = eigenvalue(i);
            if (eigenvalue(i) < 0)
                eigenStruct.t_12(i) = log(0.5)/eigenvalue(i);
            else
                eigenStruct.T2(i) = log(2)/eigenvalue(i);
            end
            eigenStruct.Tau(i) = -1/eigenvalue(i);

        else
            eigenStruct.mode{i}  = 'Oscillatory';
            eigenStruct.value(i) = eigenvalue(i);
            eigenStruct.Period(i) = 2*pi/abs(imag(eigenvalue(i)));
            if (real(eigenvalue(i)) < 0)
                eigenStruct.t_12(i) = log(0.5)/real(eigenvalue(i));
            else
                eigenStruct.T2(i) = log(2)/real(eigenvalue(i));
            end
            eigenStruct.freqNat(i) = abs(eigenvalue(i));
            eigenStruct.Damp(i) = -real(eigenvalue(i))/abs(eigenvalue(i));
        end

    end

end

