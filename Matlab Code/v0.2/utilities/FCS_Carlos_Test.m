%% CLEAR VARIABLES
clc
close all
clearvars -except LD


%% FLIGHT CONDITION (FC)
FC.hs  =   12192;                    	  % Altura en [m]
FC.a0  =    295;                          % Velocidad del sonido [m/s]
FC.us  =   206.35;                       % Velocidad condicion de referencia [m/s]
FC.Ms  =   FC.us/FC.a0;      			  % Número de Mach
FC.rho =   0.3016;					      % Densidad del aire en [kg/m^3]
FC.qs  =   0.5*FC.rho*(FC.us^2);          % Presion dinámica condicion de referencia [Pa]
FC.gravity = 9.81;                        % Gravedad en [m/s^2]
 

%% GEOMETRIC DATA (GEO)
Geo.Sw  =  21.3677;          % Superficie alar [m^2]
Geo.c   =  2.1336;          % Cuerda [m]
Geo.b   =  10.3632;         % envergadura [m]
Geo.alphaT = 0;               % alpha de los motores (tambien llamado Epsilon)
Geo.xT = 0; %Distancia de los motores al CdG proyectada en el eje Xs [m] (En este caso no es 0, pero no afecta por ser alphaT=0)
Geo.zT = 0; %Distancia de los motores al CdG proyectada en el eje Zs [m]
Geo.dT = cos(Geo.alphaT)*(Geo.zT/Geo.c) - sin(Geo.alphaT)*(Geo.xT/Geo.c); %Brazo de momentos adimensional de los motores [-]


%% WEIGHT AND BALANCE (WB)
WB.m     = 5896.7008;       % masa en [kg]
WB.Ixx_b = 37962.9026;        % Ixx en ejes cuerpo [kg*m^2]
WB.Iyy_b = 25489.37746;        % Ixx en ejes cuerpo [kg*m^2]
WB.Izz_b = 63723.44365;        % Ixx en ejes cuerpo [kg*m^2]
WB.Jxz_b = 1762.563335;        % Ixx en ejes cuerpo [kg*m^2]


%% REFERENCE CONDITION (RC)
RC.Thetas =  0;      %Theta en la condición de referencia
RC.CLs    =  0.410; %Coeficiente de sustentacion referencia
RC.CDs    =  0.0335;  %Coeficiente de resistencia referencia obtenido de la polar
RC.CTs    =  0.0335;  %Coeficiente de traccion referencia
RC.Cmas   =  0.000;  %Coeficiente de momento aerodinamico de referencia
RC.CmTs   =  0.000;  %Coeficiente de momentos debido a los motores de referencia
RC.alphas =  deg2rad(2.7);      %Ángulo entre la Xbody y Xstability en [º] es cero por definicion


%% ADIMENSIONAL STABILITY DERIVATIVES (ASD) --> LONGITUDINAL
% 0
    ASD.long.CD_0       =   0.0216;
    ASD.long.CL_0       =   0.130;
    ASD.long.Cma_0      =   0.050;
% u(Solo se tienen en cuenta si volamos en compresible)
    ASD.long.CD_u       =   0.104;
    ASD.long.CL_u       =   0.400;
    ASD.long.Cma_u      =   0.050;
% alpha
    ASD.long.CD_alpha   =   0.300;
    ASD.long.CL_alpha   =   5.840;
    ASD.long.Cma_alpha  =  -0.640;
% alphaDot
    ASD.long.CD_alphaDot  =  0.000;
    ASD.long.CL_alphaDot  =  2.200;  % CL respecto a alpha punto adimensional
    ASD.long.Cma_alphaDot = -6.700;
% q
    ASD.long.CD_q       =   0.000;
    ASD.long.CL_q       =   4.700;
    ASD.long.Cma_q      = -15.500;
% deltae
    ASD.long.CD_deltae  =  0.000;
    ASD.long.CL_deltae  =  0.460;
    ASD.long.Cma_deltae = -1.240;
% deltae_dot   
    ASD.long.CD_deltaeDot  = 0.000;
    ASD.long.CL_deltaeDot  = 0.000;
    ASD.long.Cma_deltaeDot = 0.000;   
% Motores
    ASD.long.CT_u       =  -0.07; %To be calculated in detail...
    ASD.long.CT_alpha   =   0; %To be calculated in detail...
    ASD.long.CT_deltae  =   0; %To be calculated in detail...
    ASD.long.CmT_u      =  -0.003;
    ASD.long.CmT_alpha  =   0.000;

    
%Derivadas de Estabilidad ASD longitudinales en ejes estabilidad a partir de las parciales
%(Hipotesis de alpha pequeña)
% Static
    ASD.long.Cxs      =  RC.CTs*cos(Geo.alphaT) + RC.CLs*RC.alphas - RC.CDs;
    ASD.long.Czs      = -RC.CTs*sin(Geo.alphaT) - RC.CLs - RC.CDs*RC.alphas;
    ASD.long.Cms      =  RC.Cmas + RC.CmTs;
% u
    ASD.long.Cx_u      =  ASD.long.CT_u * cos(Geo.alphaT) - ASD.long.CD_u;
    ASD.long.Cz_u      = -ASD.long.CL_u - ASD.long.CT_u * sin(Geo.alphaT);
    ASD.long.Cm_u      =  ASD.long.Cma_u + ASD.long.CmT_u;
% alpha
    ASD.long.Cx_alpha  =  ASD.long.CT_alpha*cos(Geo.alphaT) + RC.CLs - ASD.long.CD_alpha + ASD.long.CL_alpha*RC.alphas;
    ASD.long.Cz_alpha  = -ASD.long.CT_alpha*sin(Geo.alphaT) - ASD.long.CL_alpha - RC.CDs - ASD.long.CD_alpha*RC.alphas;
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
    ASD.long.Cx_deltae =  ASD.long.CT_deltae*cos(Geo.alphaT) + ASD.long.CL_deltae*RC.alphas - ASD.long.CD_deltae;
    ASD.long.Cz_deltae = -ASD.long.CT_deltae*sin(Geo.alphaT) - ASD.long.CL_deltae - ASD.long.CD_deltae*RC.alphas;
    ASD.long.Cm_deltae =  ASD.long.Cma_deltae;     
% deltae_dot   
    ASD.long.Cx_deltaeDot = -ASD.long.CD_deltaeDot;
    ASD.long.Cz_deltaeDot = -ASD.long.CL_deltaeDot;
    ASD.long.Cm_deltaeDot =  ASD.long.Cma_deltaeDot; 
    


%% ADIMENSIONAL STABILITY DERIVATIVES (ASD) --> LATERAL-DIRECTIONAL
%beta
    ASD.lat.CY_beta   =  -0.730;
    ASD.lat.Cl_beta   =  -0.110;
    ASD.lat.Cn_beta   =   0.127;
%p
    ASD.lat.CY_p      =  -0.000;
    ASD.lat.Cl_p      =  -0.450;
    ASD.lat.Cn_p      =  -0.008;
%r
    ASD.lat.CY_r      =   0.400;
    ASD.lat.Cl_r      =   0.160;
    ASD.lat.Cn_r      =  -0.200;
%deltar
    ASD.lat.CY_deltar =  -0.140;
    ASD.lat.Cl_deltar =  -0.019;
    ASD.lat.Cn_deltar =   0.074;
%deltaa
    ASD.lat.CY_deltaa =   0.000;
    ASD.lat.Cl_deltaa =   0.178;
    ASD.lat.Cn_deltaa =  -0.020;
%deltarDot
    ASD.lat.CY_deltarDot = 0.000;
    ASD.lat.Cl_deltarDot = 0.000;
    ASD.lat.Cn_deltarDot = 0.000;
%deltaaDot
    ASD.lat.CY_deltaaDot = 0.000;
    ASD.lat.Cl_deltaaDot = 0.000;
    ASD.lat.Cn_deltaaDot = 0.000;   
%Motores
    ASD.lat.CT_beta =   0;

    
%% ADIMENSIONAL STABILITY DERIVATIVES (ASD) --> PARAMETERS
   ASD.long.mu = WB.m / (0.5 * FC.rho * Geo.Sw * Geo.c);
   ASD.lat.mu  = WB.m / (0.5 * FC.rho * Geo.Sw * Geo.b);
   ASD.long.Iy = WB.Iyy_b / (FC.rho * Geo.Sw * (Geo.c/2)^3);
   ASD.lat.Ix  = WB.Ixx_b / (FC.rho * Geo.Sw * (Geo.b/2)^3);
   ASD.lat.Iz  = WB.Izz_b / (FC.rho * Geo.Sw * (Geo.b/2)^3);
   ASD.lat.Jxz = WB.Jxz_b / (FC.rho * Geo.Sw * (Geo.b/2)^3);

  
   
%% CHARACTERISTIC EQUATION (ChEq) --> LONGITUDINAL
% (Cuartica de estabilidad longitudinal)
% Metodo adimensional Longitudinal
% Construimos la matriz con las ecuaciones adimensionales hecha su transformada de Laplace
% La primera fila corresponde a "u", la segunda a "alpha" y la tercera a "theta"
% Matriz --> M = MSS*s^2 + MS*s + ML
    syms s

% MSS
    ChEq.long.MSS = zeros (3);
    ChEq.long.MSS(3,3) = -ASD.long.Iy;
    
% MS
    ChEq.long.MS  = zeros (3);
    ChEq.long.MS(1,1) =   2*ASD.long.mu;
    ChEq.long.MS(2,2) = -(2*ASD.long.mu - ASD.long.Cz_alphaDot);
    ChEq.long.MS(2,3) =   2*ASD.long.mu + ASD.long.Cz_q;
    ChEq.long.MS(3,2) =   ASD.long.Cm_alphaDot;
    ChEq.long.MS(3,3) =   ASD.long.Cm_q;

% ML
    ChEq.long.ML  = zeros (3);
    ChEq.long.ML(1,1) =  2*ASD.long.Czs*tan(RC.Thetas) - ASD.long.Cx_u;
    ChEq.long.ML(1,2) = -ASD.long.Cx_alpha;
    ChEq.long.ML(1,3) = -ASD.long.Czs;
    ChEq.long.ML(2,1) =  2*ASD.long.Czs + ASD.long.Cz_u;
    ChEq.long.ML(2,2) =  ASD.long.Cz_alpha;
    ChEq.long.ML(2,3) =  ASD.long.Czs*tan(RC.Thetas);
    ChEq.long.ML(3,1) =  ASD.long.Cm_u;
    ChEq.long.ML(3,2) =  ASD.long.Cm_alpha;
    ChEq.long.ML(3,3) =  0;

% Matriz --> M = MSS*s^2 + MS*s + ML
% Se obtiene la matriz dimensional en el dominio del tiempo sustituyendo
% por D por t*·d/dt, donde t*=c/(2us) y se dimensionaliza la velocidad
% u = us·u_adim
    ChEq.long.M = ChEq.long.MSS.*(Geo.c/(2*FC.us))^2.*s^2 + ChEq.long.MS.*(Geo.c/(2*FC.us)).*s + ChEq.long.ML;
    ChEq.long.M(1,:) = ChEq.long.M(1,:) .* FC.us;
    
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
            Modes.Phugoid.mode  = ChEq.long.EigenStruct.mode{I(1)};
            Modes.Phugoid.value = [ChEq.long.EigenStruct.value(I(1)),ChEq.long.EigenStruct.value(I(2))];
            Modes.Phugoid.freqNat = ChEq.long.EigenStruct.freqNat(I(1));
            Modes.Phugoid.Damp = ChEq.long.EigenStruct.Damp(I(1));
            Modes.Phugoid.Period = ChEq.long.EigenStruct.Period(I(1));
            Modes.Phugoid.t_12 = ChEq.long.EigenStruct.t_12(I(1));
            Modes.Phugoid.T2 = ChEq.long.EigenStruct.T2(I(1));
            Modes.Phugoid.Tau = ChEq.long.EigenStruct.Tau(I(1));

            %Short Period --> Corto Periodo
            Modes.ShortPeriod.mode  = ChEq.long.EigenStruct.mode{I(3)};
            Modes.ShortPeriod.value = [ChEq.long.EigenStruct.value(I(3)),ChEq.long.EigenStruct.value(I(4))];
            Modes.ShortPeriod.freqNat = ChEq.long.EigenStruct.freqNat(I(3));
            Modes.ShortPeriod.Damp = ChEq.long.EigenStruct.Damp(I(3));
            Modes.ShortPeriod.Period = ChEq.long.EigenStruct.Period(I(3));
            Modes.ShortPeriod.t_12 = ChEq.long.EigenStruct.t_12(I(3));
            Modes.ShortPeriod.T2 = ChEq.long.EigenStruct.T2(I(3));
            Modes.ShortPeriod.Tau = ChEq.long.EigenStruct.Tau(I(3));
            clear I
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
% La primera fila corresponde a "beta", la segunda a "Phi" y la tercera a "r"
% Matriz --> M = MSS*s^2 + MS*s + ML
    syms s

% MSS
    ChEq.lat.MSS = zeros (3);
    ChEq.lat.MSS(2,2) =  ASD.lat.Ix;
    ChEq.lat.MSS(3,2) = -ASD.lat.Jxz;
    
% MS
    ChEq.lat.MS  = zeros (3);
    ChEq.lat.MS(1,1) =   2*ASD.lat.mu;
    ChEq.lat.MS(1,2) = - ASD.lat.CY_p;
    ChEq.lat.MS(2,2) = - ASD.lat.Cl_p;
    ChEq.lat.MS(2,3) = - ASD.lat.Jxz;
    ChEq.lat.MS(3,2) = - ASD.lat.Cn_p;
    ChEq.lat.MS(3,3) =   ASD.lat.Iz;

% ML
    ChEq.lat.ML  = zeros (3);
    ChEq.lat.ML(1,1) = - ASD.lat.CY_beta;
    ChEq.lat.ML(1,2) = - RC.CLs;
    ChEq.lat.ML(1,3) =   2*ASD.lat.mu - ASD.lat.CY_r;
    ChEq.lat.ML(2,1) = - ASD.lat.Cl_beta;
    ChEq.lat.ML(2,2) =   0;
    ChEq.lat.ML(2,3) = - ASD.lat.Cl_r;
    ChEq.lat.ML(3,1) = - ASD.lat.Cn_beta;
    ChEq.lat.ML(3,2) =   0;
    ChEq.lat.ML(3,3) = - ASD.lat.Cn_r;

% Matriz --> M = MSS*s^2 + MS*s + ML
% Se obtiene la matriz dimensional en el dominio del tiempo sustituyendo
% por D por t*·d/dt, donde t*=b/(2us) y se dimensionaliza la velocidad
% u = us·u_adim
    ChEq.lat.M = ChEq.lat.MSS.*(Geo.b/(2*FC.us))^2.*s^2 + ChEq.lat.MS.*(Geo.b/(2*FC.us)).*s + ChEq.lat.ML;
    ChEq.lat.M(3,:) = ChEq.lat.M(3,:) .* (2*FC.us/Geo.b);

% Calculamos el determinante y los coeficientes de este
% Segun el criterio de Routh, si la cuartica de estabilidad es:
% A*s^4 + B*s^3 + C*s^2 + D*s + E = 0
% El sistema es estable cuando: A>0, B>0, C>0, D>0, E>0 y [D(BC-AD)-(B^2)E]>0
    ChEq.lat.det_M  = det(ChEq.lat.M);
    ChEq.lat.Coeffs = sym2poly(ChEq.lat.det_M); 
    
% Calculamos los autovalores e identificamos el modo fugoide y corto periodo
    ChEq.lat.Eigenvalues = roots(ChEq.lat.Coeffs);
    if (length(ChEq.lat.Eigenvalues) == 4) 
        for j=1:4; R(j)=isreal(ChEq.lat.Eigenvalues(j)); end; clear j %#ok<SAGROW>
        if (sum(R) == 2)
        	ChEq.lat.EigenStruct = getEigenData(ChEq.lat.Eigenvalues);
            [~,I] = sort(abs(real(ChEq.lat.EigenStruct.value)));
            %Spiral --> Modo Espiral
            Modes.Spiral.mode  = ChEq.lat.EigenStruct.mode{I(1)};
            Modes.Spiral.value = ChEq.lat.EigenStruct.value(I(1));
            Modes.Spiral.freqNat = ChEq.lat.EigenStruct.freqNat(I(1));
            Modes.Spiral.Damp = ChEq.lat.EigenStruct.Damp(I(1));
            Modes.Spiral.Period = ChEq.lat.EigenStruct.Period(I(1));
            Modes.Spiral.t_12 = ChEq.lat.EigenStruct.t_12(I(1));
            Modes.Spiral.T2 = ChEq.lat.EigenStruct.T2(I(1));
            Modes.Spiral.Tau = ChEq.lat.EigenStruct.Tau(I(1));
            
            %Roll --> Modo de Convergencia en Balance
            Modes.Roll.mode  = ChEq.lat.EigenStruct.mode{I(4)};
            Modes.Roll.value = ChEq.lat.EigenStruct.value(I(4));
            Modes.Roll.freqNat = ChEq.lat.EigenStruct.freqNat(I(4));
            Modes.Roll.Damp = ChEq.lat.EigenStruct.Damp(I(1));
            Modes.Roll.Period = ChEq.lat.EigenStruct.Period(I(4));
            Modes.Roll.t_12 = ChEq.lat.EigenStruct.t_12(I(4));
            Modes.Roll.T2 = ChEq.lat.EigenStruct.T2(I(4));
            Modes.Roll.Tau = ChEq.lat.EigenStruct.Tau(I(4));
            
            %Dutch Roll --> Modo de Balanceo Holandes
            Modes.DutchRoll.mode  = ChEq.lat.EigenStruct.mode{I(2)};
            Modes.DutchRoll.value = [ChEq.lat.EigenStruct.value(I(2)),ChEq.long.EigenStruct.value(I(3))];
            Modes.DutchRoll.freqNat = ChEq.lat.EigenStruct.freqNat(I(2));
            Modes.DutchRoll.Damp = ChEq.lat.EigenStruct.Damp(I(2));
            Modes.DutchRoll.Period = ChEq.lat.EigenStruct.Period(I(2));
            Modes.DutchRoll.t_12 = ChEq.lat.EigenStruct.t_12(I(2));
            Modes.DutchRoll.T2 = ChEq.lat.EigenStruct.T2(I(2));
            Modes.DutchRoll.Tau = ChEq.lat.EigenStruct.Tau(I(2));
            clear I
        else
            wrn = msgbox({'Los modos lateral-direccionales no son 2 reales y 2 osscilatorios.',...
                          'Se recomienda revisar a mano lo que está sucediendo.'},...
                          'Warning','warn');
            uiwait(wrn);
            disp('Los modos lateral-direccionales no son 2 reales y 2 osscilatorios.')
            disp('Se recomienda revisar a mano lo que está sucediendo.')
            pause
        end
        clear R
    else
        wrn = msgbox({'El número de modos lateral-direccionales es diferente de 4.',...
                      'Se recomienda revisar a mano lo que está sucediendo.'},...
                      'Warning','warn');
        uiwait(wrn);
        disp('El número de modos lateral-direccionales es diferente de 4.')
        disp('Se recomienda revisar a mano lo que está sucediendo.')
        pause
    end
    
    
    
%% TRANSFER FUNCTIONS (TF) --> LONGITUDINAL

% FUNCIÓN DE TRANSFERENCIA DE ELEVADOR A VELOCIDAD

% Construimos el numerador
% MK2
    TF.MK2 = ChEq.long.MSS;
    TF.MK2(:,1) = 0 ;
% MK1
    TF.MK1 = ChEq.long.MS;
    TF.MK1(:,1)  = 0 ;
    TF.MK1(3,1)  = ASD.long.Cm_deltaeDot;
% MK0
    TF.MK0 = ML;
    TF.MK0(:,1)   = 0 ;
    TF.MK0(1,1)   = ASD.long.Cx_deltae;
    TF.MK0(2,1)   = ASD.long.Cz_deltae;
    TF.MK0(3,1)   = ASD.long.Cm_deltae;
% Matriz --> MK = MK2*s^2 + MK1*s + MK0
% Se obtiene la matriz dimensional en el dominio del tiempo sustituyendo
% por D por t*·d/dt, donde t*=c/(2us) y se dimensionaliza la velocidad
% u = us·u_adim
    TF.MK = TF.MK2.*(Geo.c/(2*FC.us))^2.*s^2 + TF.MK1.*(Geo.c/(2*FC.us)).*s + TF.MK0;
    TF.MK(3,:) = TF.MK(3,:) .* (2*FC.us/Geo.c);



roskam.Elev_vel.sist    = simplify(s^2*KKK+s*KK+K);
roskam.Elev_vel.ec      = det(roskam.Elev_vel.sist);
roskam.Elev_vel.coef    = sym2poly(roskam.Elev_vel.ec);
roskam.Elev_vel.raices  = roots(roskam.Elev_vel.coef);
roskam.Elev_vel.tau1    = -1/roskam.Elev_vel.raices(3);
roskam.Elev_vel.tau2    = -1/roskam.Elev_vel.raices(2);
roskam.Elev_vel.K       =  roskam.Elev_vel.coef(1)*(-1/roskam.Elev_vel.tau1)*(-1/roskam.Elev_vel.tau2)/((-1/roskam.cuart.tau1_SP)*(-1/roskam.cuart.tau2_SP)*roskam.cuart.w_nF^2*roskam.cuart.Cuarticalong_coef(1));

% Elevador a Ángulo de Ataque
K         = ML;
KK        = MS;
KKK       = MSS;
K(:,2)    = 0 ;
KK(:,2)   = 0 ;
KKK(:,2)  = 0 ;
KK(1,2)   = roskam.X_deltae;
KK(2,2)   = roskam.Z_deltae;
KK(3,2)   = roskam.M_deltae;

roskam.Elev_ataque.sist    = simplify(s^2*KKK+s*KK+K);
roskam.Elev_ataque.ec      = det(roskam.Elev_ataque.sist);
roskam.Elev_ataque.coef    = sym2poly(roskam.Elev_ataque.ec);
roskam.Elev_ataque.raices  = roots(roskam.Elev_ataque.coef);
roskam.Elev_ataque.Wn      = abs(roskam.Elev_ataque.raices(3));
roskam.Elev_ataque.Tsi     =-real(roskam.Elev_ataque.raices(3))/roskam.Elev_ataque.Wn;
roskam.Elev_ataque.Tau     = -1/roskam.Elev_ataque.raices(2);
roskam.Elev_ataque.K       = -roskam.Elev_ataque.coef(1)*roskam.Elev_ataque.Wn^2*(-1/roskam.Elev_ataque.Tau)/((-1/roskam.cuart.tau1_SP)*(-1/roskam.cuart.tau2_SP)*roskam.cuart.w_nF^2*roskam.cuart.Cuarticalong_coef(1));

% Elevador a Ángulo de Asiento
K         = ML;
KK        = MS;
KKK       = MSS;
K(:,3)    = 0 ;
KK(:,3)   = 0 ;
KKK(:,3)  = 0 ;
KK(1,3)   = roskam.X_deltae;
KK(2,3)   = roskam.Z_deltae;
KK(3,3)   = roskam.M_deltae;

roskam.Elev_asiento.sist    = simplify(s^2*KKK+s*KK+K);
roskam.Elev_asiento.ec      = det(roskam.Elev_asiento.sist);
roskam.Elev_asiento.coef    = sym2poly(roskam.Elev_asiento.ec);
roskam.Elev_asiento.raices  = roots(roskam.Elev_asiento.coef);
roskam.Elev_asiento.Tau1     = -1/roskam.Elev_asiento.raices(2);
roskam.Elev_asiento.Tau2     = -1/roskam.Elev_asiento.raices(3);
roskam.Elev_asiento.K        = roskam.Elev_asiento.coef(1)*(-1/roskam.Elev_asiento.Tau1)*(-1/roskam.Elev_asiento.Tau2)/((-1/roskam.cuart.tau1_SP)*(-1/roskam.cuart.tau2_SP)*roskam.cuart.w_nF^2*roskam.cuart.Cuarticalong_coef(1));
