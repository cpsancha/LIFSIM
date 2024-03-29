%% CLEAR VARIABLES
clc
close all
clearvars -except LD

%% CONVERSI�N DE UNIDADES A S.I.
% Definicion de los parametros de la practica: LIBIS Cruise configuration

%Definimos los cambios de unidades
% ft2m   =   1/3.28084;      %Factor de conversion pies a metros 
% lb2kg  =   0.453592;       %Factor de conversion libra a kg
% lbft2Pa =  47.88;          %Factor de conversion presion din�mica a pascales
% slgft2kgm = 1.355817962;   %Factor de conversion slug*pie^2 a kg*m^2
% deg2rad   = pi/180;        %Factor de conversion grados a radianes

%% Flight Condition (FC)
FC.hs  =   91.44;                    	  % Altura en [m]
FC.a0  =    340;                          % Velocidad del sonido [m/s]
FC.us  =   42.22;                         % Velocidad condicion de referencia [m/s]
FC.Ms  =   FC.us/FC.a0;      			  % N�mero de Mach
FC.rho =   1.21;					      % Densidad del aire en [kg/m^3]
FC.qs  =   0.5*FC.rho*(FC.us^2);          % Presion din�mica condicion de referencia [Pa]
FC.gravity   = 9.81;                      % Gravedad en [m/s^2]
FC.alfa_xbxs = 0;                         % �ngulo entre la Xbody y Xstability en [�]
 
%% Datos Geometricos (GEO)
Geo.Sw  =  LD.sref;          % Superficie alar [m^2]
Geo.c   =  LD.cbar;          % Cuerda [m]
Geo.b   =  LD.blref;         % envergadura [m]
Geo.alfaT = 0;               % alfa de los motores

%% Pesos y Equilibrado (WB) (weight and balance)
WB.m     = LD.Inertia.mass;       % masa en [kg]
WB.Ixx_b = LD.Inertia.Ixx;        % Ixx en ejes cuerpo [kg*m^2]
WB.Iyy_b = LD.Inertia.Iyy;        % Ixx en ejes cuerpo [kg*m^2]
WB.Izz_b = LD.Inertia.Izz;        % Ixx en ejes cuerpo [kg*m^2]
WB.Ixz_b = LD.Inertia.Ixz;        % Ixx en ejes cuerpo [kg*m^2]

%% Condicion de Referencia RC
RC.Thetas = 0;      %Theta en la condici�n de referencia
RC.Cls    = WB.m*FC.gravity*cos(RC.Thetas)/(FC.qs*Geo.Sw);   %Coeficiente de sustentacion referencia
RC.Cds    = 0.0262; %Coeficiente de resistencia referencia
RC.Ct_xs  = 0.0262; %Coeficiente de traccion referencia en eje Xs
RC.Cms    = 0;      %Coeficiente de momento referencia
RC.Cmts   = 0;      %Coeficiente de momento de traccion referencia

%% Derivadas de Estabilidad Adimensionales longitudimales ASD (Stability Derivatives)
% 0
    ASD.long.CD_0      =   0.024;
    ASD.long.CL_0      =   0.181;
    ASD.long.Cm_0      =   -0.003;
% u(Solo se tienen en cuenta si volamos en compresible)
    ASD.long.CD_u      =   0;
    ASD.long.CL_u      =   0;
    ASD.long.Cm_u      =   0;
% alfa
    ASD.long.CD_alfa   =   0.239;
    ASD.long.CL_alfa   =   4.138;
    ASD.long.Cm_alfa   =   -0.707;
% alfap
    ASD.long.CD_alfap  =   0;
    ASD.long.CL_alfap  =   7.411;  % CL respecto a alfa punto adimensional
    ASD.long.Cm_alfap  =  -4.063;
% q
    ASD.long.CD_q      =   0;
    ASD.long.CL_q      =   6.248;
    ASD.long.Cm_q      =  -18.410;
% deltae
    ASD.long.CD_deltae =   0.005;
    ASD.long.CL_deltae =   1.391;
    ASD.long.Cm_deltae =  -2.930;
% Motores
    ASD.long.CT_u     =  0; %To be calculated in detail...
    ASD.long.Cm_Tu    =  0; %Los motores no producen momento al acelerar
    ASD.long.CT_alfa  =  0; %To be calculated in detail...
    ASD.long.Cm_Talfa =  0; %Los motores no producen momento al ser 0 el brazo

    
%% Derivadas de Estabilidad Adimensionales lateral-direccionales ASD (Stability Derivatives)
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
%Motores
    ASD.lat.Cn_T_beta =   0;

    
    
%% DERIVADAS DE ESTABILIDAD DIMENSIONALES
% Una vez se tienen los datos, se pasa a obtener las derivadas de estabilidad dimensionales

% [T, a, P, FC.rho] = atmosisa(FC.hs);

%Adimensionalizadores longitudinales
Adimlong_1 = FC.rho*FC.us*Geo.Sw;             %(rho*Us*S)
Adimlong_2 = FC.rho*FC.us*Geo.Sw/2;           %(rho*Us*S)/2
Adimlong_3 = FC.rho*(FC.us)^2*Geo.Sw/2;       %(rho*(Us)^2*S)/2
Adimlong_4 = FC.rho*Geo.Sw*Geo.c/4;           %(rho*S*c)/4
Adimlong_5 = FC.rho*FC.us*Geo.Sw*Geo.c/4;     %(rho*Us*S*c)/4
Adimlong_6 = FC.rho*FC.us*Geo.Sw*Geo.c/2;     %(rho*Us*S*c)/2
Adimlong_7 = FC.rho*Geo.Sw*(Geo.c)^2/4;       %(rho*S*c^2)/4
Adimlong_8 = FC.rho*FC.us*Geo.Sw*(Geo.c)^2/4; %(rho*Us*S*c^2)/4
Adimlong_9 = FC.rho*(FC.us)^2*Geo.Sw*Geo.c/2; %(rho*(Us)^2*S*c)/2
Adimlong_10= FC.rho*(FC.us)^2*Geo.Sw*Geo.c/4; %(rho*(Us)^2*S*c)/4

%Derivadas de Estabilidad ASD longitudinales a partir de las que se tienen
%Sin Hipotesis de alfa peque�a 
% Static
ASD.long.Coeffs.C_xs     =    RC.Ct_xs*cos(Geo.alfaT)+RC.Cls*FC.alfa_xbxs-RC.Cds;
ASD.long.Coeffs.C_zs     =  - RC.Cls - RC.Cds*FC.alfa_xbxs - RC.Ct_xs*sin(Geo.alfaT);  %RC.Ct_xs
ASD.long.Coeffs.C_ms     =    RC.Cms - RC.Cmts;
% u
ASD.long.Coeffs.C_xu     =  ASD.long.CT_u*cos(Geo.alfaT) + ASD.long.CL_u*FC.alfa_xbxs- ASD.long.CD_u;
ASD.long.Coeffs.C_zu     = -ASD.long.CL_u - ASD.long.CD_u*FC.alfa_xbxs - ASD.long.CT_u*sin(Geo.alfaT);
ASD.long.Coeffs.C_mu     =  ASD.long.Cm_u- ASD.long.Cm_Tu ;
% alfa
ASD.long.Coeffs.C_xalfa  =  ASD.long.CT_alfa*cos(Geo.alfaT) + RC.Cls - ASD.long.CD_alfa + ASD.long.CL_alfa*FC.alfa_xbxs;
ASD.long.Coeffs.C_zalfa  = -ASD.long.CT_alfa*sin(Geo.alfaT) - ASD.long.CL_alfa - RC.Cds - ASD.long.CD_alfa*FC.alfa_xbxs;
ASD.long.Coeffs.C_malfa  =  ASD.long.Cm_alfa - ASD.long.Cm_Talfa;
% alfap
ASD.long.Coeffs.C_xalfap = -ASD.long.CD_alfap;
ASD.long.Coeffs.C_zalfap = -ASD.long.CL_alfap;
ASD.long.Coeffs.C_malfap =  ASD.long.Cm_alfap;
% q
ASD.long.Coeffs.C_xq     = -ASD.long.CD_q;
ASD.long.Coeffs.C_zq     = -ASD.long.CL_q;
ASD.long.Coeffs.C_mq     =  ASD.long.Cm_q;
% deltae
ASD.long.Coeffs.C_xdeltae = -ASD.long.CD_deltae;
ASD.long.Coeffs.C_zdeltae = -ASD.long.CL_deltae;
ASD.long.Coeffs.C_mdeltae =  ASD.long.Cm_deltae;
  
%Dimensionalizar las derivadas
SD.long.Xu      = Adimlong_1*ASD.long.Coeffs.C_xs   + Adimlong_2*ASD.long.Coeffs.C_xu; %X_u 
SD.long.Zu      = Adimlong_1*ASD.long.Coeffs.C_zs   + Adimlong_2*ASD.long.Coeffs.C_zu; %Z_u
SD.long.Mu      = Adimlong_6*2*ASD.long.Coeffs.C_ms + Adimlong_6*ASD.long.Coeffs.C_mu; %M_u
SD.long.Xw      = Adimlong_2*ASD.long.Coeffs.C_xalfa;                  %X_w
SD.long.Zw      = Adimlong_2*ASD.long.Coeffs.C_zalfa;                  %Z_w
SD.long.Mw      = Adimlong_6*ASD.long.Coeffs.C_malfa;                  %M_w
SD.long.Zwp     = Adimlong_4*ASD.long.Coeffs.C_zalfap;                 %Z_wp
SD.long.Mwp     = Adimlong_7*ASD.long.Coeffs.C_malfap;                 %M_wp
SD.long.Zq      = Adimlong_5*ASD.long.Coeffs.C_zq;                     %Z_q
SD.long.Mq      = Adimlong_8*ASD.long.Coeffs.C_mq;                     %M_q
SD.long.Xdeltae = Adimlong_3*ASD.long.Coeffs.C_xdeltae;                %X_deltae
SD.long.Zdeltae = Adimlong_3*ASD.long.Coeffs.C_zdeltae;                %Z_deltae
SD.long.Mdeltae = Adimlong_9*ASD.long.Coeffs.C_mdeltae;                %M_deltae


%Adimensionalizadores lateral-direccionales
Adimlat_1 = FC.rho*FC.us*Geo.Sw/2;               %(FC.rho*Us*S)/2
Adimlat_2 = FC.rho*(FC.us)*Geo.Sw*Geo.b/2;       %(FC.rho*Us*S*b)/2
Adimlat_3 = FC.rho*(FC.us)*Geo.Sw*Geo.b/4;       %(FC.rho*Us*S*b)/4
Adimlat_4 = FC.rho*(FC.us)*Geo.Sw*(Geo.b)^2/4;   %(FC.rho*Us*S*b^2)/4
Adimlat_5 = FC.rho*(FC.us)^2*Geo.Sw*Geo.b/2;     %(FC.rho*Us^2*S*b)/2
Adimlat_6 = FC.rho*(FC.us)^2*Geo.Sw/2;           %(FC.rho*Us^2*S)/2

%Derivadas de Estabilidad SD
SD.lat.Yv      = Adimlat_1*ASD.lat.CY_beta;         %Y_v
SD.lat.Lv      = Adimlat_2*ASD.lat.Cl_beta;         %L_v
SD.lat.Nv      = Adimlat_2*ASD.lat.Cn_beta;         %N_v
SD.lat.Yp      = Adimlat_3*ASD.lat.CY_p;            %Y_p
SD.lat.Lp      = Adimlat_4*ASD.lat.Cl_p;            %L_p
SD.lat.Np      = Adimlat_4*ASD.lat.Cn_p;            %N_p
SD.lat.Yr      = Adimlat_3*ASD.lat.CY_r;            %Y_r
SD.lat.Lr      = Adimlat_4*ASD.lat.Cl_r;            %L_r
SD.lat.Nr      = Adimlat_4*ASD.lat.Cn_r;            %N_r
SD.lat.Ydeltar = Adimlat_6*ASD.lat.CY_deltar;       %Y_deltar
SD.lat.Ldeltar = Adimlat_5*ASD.lat.Cl_deltar;       %L_deltar
SD.lat.Ndeltar = Adimlat_5*ASD.lat.Cn_deltar;       %N_deltar
SD.lat.Ldeltaa = Adimlat_5*ASD.lat.Cl_deltaa;       %L_deltaa
SD.lat.Ndeltaa = Adimlat_5*ASD.lat.Cn_deltaa;       %N_deltaa

clear -regexp ^Adimlong_ ^Adimlat_



%% SISTEMA ROSKAM
roskam.X_u      = -FC.qs*Geo.Sw*(ASD.long.CD_u+2*RC.Cds)/(WB.m*FC.us);
roskam.X_Tu     =  FC.qs*Geo.Sw*(ASD.long.CT_u+2*RC.Ct_xs)/(WB.m*FC.us);
roskam.X_alfa   = -FC.qs*Geo.Sw*(ASD.long.CD_alfa-RC.Cls)/WB.m;
roskam.X_deltae = -FC.qs*Geo.Sw*ASD.long.CD_deltae/WB.m;
roskam.Z_u      = -FC.qs*Geo.Sw*(ASD.long.CL_u+2*RC.Cls)/(WB.m*FC.us);
roskam.Z_alfa   = -FC.qs*Geo.Sw*(ASD.long.CL_alfa-RC.Cds)/WB.m;
roskam.Z_alfap  = -FC.qs*Geo.Sw*Geo.c*ASD.long.CL_alfap/(2*WB.m*FC.us);
roskam.Z_q      = -FC.qs*Geo.Sw*Geo.c*ASD.long.CL_q/(2*WB.m*FC.us);
roskam.Z_deltae = -FC.qs*Geo.Sw*ASD.long.CL_deltae/WB.m;
roskam.M_u      =  FC.qs*Geo.Sw*Geo.c*(ASD.long.Cm_u+2*RC.Cms)/(WB.Iyy_b*FC.us);
roskam.M_Tu     =  FC.qs*Geo.Sw*Geo.c*(ASD.long.Cm_Tu+2*RC.Cmts)/(WB.Iyy_b*FC.us);
roskam.M_alfa   =  FC.qs*Geo.Sw*Geo.c*ASD.long.Cm_alfa/WB.Iyy_b;
roskam.M_Talfa  =  FC.qs*Geo.Sw*Geo.c*ASD.long.Cm_Talfa/WB.Iyy_b;
roskam.M_alfap  =  FC.qs*Geo.Sw*Geo.c^2*ASD.long.Cm_alfap/(2*WB.Iyy_b*FC.us);
roskam.M_q      =  FC.qs*Geo.Sw*Geo.c^2*ASD.long.Cm_q/(2*WB.Iyy_b*FC.us);
roskam.M_deltae =  FC.qs*Geo.Sw*Geo.c*ASD.long.Cm_deltae/WB.Iyy_b;
syms s


%Matriz MS
MSS = zeros (3);
MS  = zeros (3);
ML  = zeros (3);

MSS(3,3) = 1;

MS(1,1) = 1;
MS(2,2) = (FC.us-roskam.Z_alfap);
MS(2,3) = -(roskam.Z_q+FC.us);
MS(3,2) = -roskam.M_alfap;
MS(3,3) = -roskam.M_q;

ML(1,1) = -roskam.X_u-roskam.X_Tu;
ML(1,2) = -roskam.X_alfa;
ML(1,3) =  FC.gravity*cos(RC.Thetas);
ML(2,1) = -roskam.Z_u;
ML(2,2) = -roskam.Z_alfa;
ML(2,3) = FC.gravity*sin(RC.Thetas);
ML(3,1) = -(roskam.M_u+roskam.M_Tu);
ML(3,2) = -(roskam.M_alfa+roskam.M_Talfa);

roskam.cuart.matriz            = s^2*MSS+s*MS+ML;
roskam.cuart.deter_matriz      = det(roskam.cuart.matriz);
roskam.cuart.Cuarticalong_coef = sym2poly(roskam.cuart.deter_matriz); 
roskam.cuart.Autovalores       = roots(roskam.cuart.Cuarticalong_coef);
roskam.cuart.tau1_SP           = -1/roskam.cuart.Autovalores(1); 
roskam.cuart.tau2_SP           = -1/roskam.cuart.Autovalores(2); 
roskam.cuart.w_nF              = abs(roskam.cuart.Autovalores(3));
roskam.cuart.tsiF               =-real(roskam.cuart.Autovalores(3))/roskam.cuart.w_nF ;


%% CRAMER de cada modo

% Elevador a velocidad
K         = ML;
KK        = MS;
KKK       = MSS;
K(:,1)    = 0 ;
KK(:,1)   = 0 ;
KKK(:,1)  = 0 ;
KK(1,1)   = roskam.X_deltae;
KK(2,1)   = roskam.Z_deltae;
KK(3,1)   = roskam.M_deltae;

roskam.Elev_vel.sist    = simplify(s^2*KKK+s*KK+K);
roskam.Elev_vel.ec      = det(roskam.Elev_vel.sist);
roskam.Elev_vel.coef    = sym2poly(roskam.Elev_vel.ec);
roskam.Elev_vel.raices  = roots(roskam.Elev_vel.coef);
roskam.Elev_vel.tau1    = -1/roskam.Elev_vel.raices(3);
roskam.Elev_vel.tau2    = -1/roskam.Elev_vel.raices(2);
roskam.Elev_vel.K       =  roskam.Elev_vel.coef(1)*(-1/roskam.Elev_vel.tau1)*(-1/roskam.Elev_vel.tau2)/((-1/roskam.cuart.tau1_SP)*(-1/roskam.cuart.tau2_SP)*roskam.cuart.w_nF^2*roskam.cuart.Cuarticalong_coef(1));

% Elevador a �ngulo de Ataque
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

% Elevador a �ngulo de Asiento
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


%% Metodo adimensional Longitudinal
%Construimos matriz con las ecuaciones adimensionales dimensionalizada y hecha su transformada de Laplace
syms s

%Matriz MS
MSS = zeros (4);
MS = zeros (4);
ML = zeros (4);
MS(1,1) =  (2*WB.m/(FC.rho*Geo.Sw*FC.us^2));
ML(1,1) = -(ASD.long.Coeffs.C_xu/FC.us);
ML(1,2) = -ASD.long.Coeffs.C_xalfa;%%
ML(1,3) = -ASD.long.Coeffs.C_zs;
ML(2,1) = -((ASD.long.Coeffs.C_zu+2*ASD.long.Coeffs.C_zs)/FC.us);
MS(2,2) =  (((4*WB.m/(FC.rho*Geo.Sw*Geo.c))-ASD.long.Coeffs.C_zalfap)*Geo.c)/(2*FC.us);%%
ML(2,2) = -ASD.long.Coeffs.C_zalfa;%%
MS(2,3) = -(((4*WB.m/(FC.rho*Geo.Sw*Geo.c))+ASD.long.Coeffs.C_zq)*Geo.c)/(2*FC.us);
ML(3,1) = -ASD.long.Coeffs.C_mu/FC.us;
ML(3,2) = -ASD.long.Coeffs.C_malfa;%%
MS(3,2) = -ASD.long.Coeffs.C_malfap*Geo.c/(2*FC.us);%%
MS(3,3) = -ASD.long.Coeffs.C_mq*Geo.c/(2*FC.us);
MS(3,4) =  WB.Iyy_b/(FC.rho*Geo.Sw*(Geo.c/2))*(1/FC.us^2);
MSS(3,3) = 0;
MS (4,3)=  Geo.c/(2*FC.us);
ML (4,4)= -Geo.c/(2*FC.us);

Cramer.cuart.matriz            = s^2*MSS+s*MS+ML;
Cramer.cuart.deter_matriz      = det(Cramer.cuart.matriz);
Cramer.cuart.Cuarticalong_coef = sym2poly(Cramer.cuart.deter_matriz); 
Cramer.cuart.Autovalores       = roots(Cramer.cuart.Cuarticalong_coef);
Cramer.cuart.tau1_SP           = 1/(-Cramer.cuart.Autovalores(1));
Cramer.cuart.tau2_SP           = 1/(-Cramer.cuart.Autovalores(2));
Cramer.cuart.w_nF              = abs(Cramer.cuart.Autovalores(3));
Cramer.cuart.tsiF              =-real(Cramer.cuart.Autovalores(3))/Cramer.cuart.w_nF ;
% No tiene corto periodo;





%% CRAMER de cada modo

% Elevador a velocidad
K  = MS;
KK = ML;
K(:,1) =0 ;
KK(1,1)= ASD.long.Coeffs.C_xdeltae;
KK(2,1)= ASD.long.Coeffs.C_zdeltae;
KK(3,1)= ASD.long.Coeffs.C_mdeltae;

Cramer.Elev_vel.sist    = simplify(s^2*MSS+s*K+KK);
Cramer.Elev_vel.ec      = det(Cramer.Elev_vel.sist);
Cramer.Elev_vel.coef    = sym2poly(Cramer.Elev_vel.ec);
Cramer.Elev_vel.raices  = roots(Cramer.Elev_vel.coef);


% Elevador a �ngulo de Ataque
K  = MS;
KK = ML;
K(:,2) =0 ;
KK(1,2)= ASD.long.Coeffs.C_xdeltae;
KK(2,2)= ASD.long.Coeffs.C_zdeltae;
KK(3,2)= ASD.long.Coeffs.C_mdeltae;

Cramer.Elev_ataque.sist    = simplify(s^2*MSS+s*K+KK);
Cramer.Elev_ataque.ec      = det(Cramer.Elev_ataque.sist);
Cramer.Elev_ataque.coef    = sym2poly(Cramer.Elev_ataque.ec);
Cramer.Elev_ataque.raices  = roots(Cramer.Elev_ataque.coef);
Cramer.Elev_ataque.Wn      = abs(Cramer.Elev_ataque.raices(2));
Cramer.Elev_ataque.Tsi     =-real(Cramer.Elev_ataque.raices(2))/Cramer.Elev_ataque.Wn;
Cramer.Elev_ataque.Tau     = -1/Cramer.Elev_ataque.raices(1);
Cramer.Elev_ataque.K       = -Cramer.Elev_ataque.coef(1)*Cramer.Elev_ataque.Wn^2*(-1/Cramer.Elev_ataque.Tau)/((-1/Cramer.cuart.tau1_SP)*(-1/Cramer.cuart.tau2_SP)*Cramer.cuart.w_nF^2*Cramer.cuart.Cuarticalong_coef(1));


% Elevador a �ngulo de Asiento
K  = MS;
KK = ML;
KKK =MSS;
K(:,3) =0 ;
KK(1,3)= ASD.long.Coeffs.C_xdeltae;
KK(2,3)= ASD.long.Coeffs.C_zdeltae;
KK(3,3)= ASD.long.Coeffs.C_mdeltae;
KKK(:,3) = 0;

Cramer.Elev_asiento.sist    = simplify(s^2*KKK+s*K+KK);
Cramer.Elev_asiento.ec      = det(Cramer.Elev_asiento.sist);
Cramer.Elev_asiento.coef    = sym2poly(Cramer.Elev_asiento.ec);
Cramer.Elev_asiento.raices  = roots(Cramer.Elev_asiento.coef);
Cramer.Elev_asiento.Tau1     = -1/Cramer.Elev_asiento.raices(1);
Cramer.Elev_asiento.Tau2     = -1/Cramer.Elev_asiento.raices(2);
Cramer.Elev_asiento.K        = Cramer.Elev_asiento.coef(1)*(-1/Cramer.Elev_asiento.Tau1)*(-1/Cramer.Elev_asiento.Tau2)/((-1/Cramer.cuart.tau1_SP)*(-1/Cramer.cuart.tau2_SP)*Cramer.cuart.w_nF^2*Cramer.cuart.Cuarticalong_coef(1));



%% Una vez tenemos las derivada de estabilidad dimensionales. Construimos el sistema de ecuaciones
%Primero el longitudinal
%matriz E*DX=A'*X+ B'*U -----> DX = A*X +B*U
E         =  zeros(4);
E(1,1)    = -WB.m;
E(2,2)    =  SD.long.Zwp - WB.m;
E(2,3)    =  0;
E(3,2)    =  SD.long.Mwp;
E(3,3)    =  SD.long.Mq;
E(3,4)    = -WB.Iyy_b;
E(4,3)    =  1;

A_p       =  zeros(4);
A_p(1,1)  =  SD.long.Xu;
A_p(1,2)  =  SD.long.Xw;
A_p(1,3)  = -WB.m*FC.gravity*cos(RC.Thetas);
A_p(2,1)  =  SD.long.Zu;
A_p(2,2)  =  SD.long.Zw;
A_p(2,3)  = -WB.m*FC.gravity*sin(RC.Thetas);
A_p(2,4)  =  SD.long.Zq + WB.m*FC.us;
A_p(3,1)  =  SD.long.Mu;
A_p(3,2)  =  SD.long.Mw;
A_p(3,4)  =  0;
A_p(4,4)  =  -1;
A_p=-A_p;

B_p       =  zeros(4,1);
B_p(1,1)  = -SD.long.Xdeltae;
B_p(2,1)  = -SD.long.Zdeltae;
B_p(3,1)  = -SD.long.Mdeltae;


A = E\A_p; %Inv(E)*A_p
B = E\B_p; %Inv(E)*B_p
C = eye(4);
D = zeros(4,1); 


% Funciones de transferencia longitudinales
[TF.long.Nums,TF.long.Den] = ss2tf(A,B,C,D);


%% Ahora el lateral-direccional
%matriz E_lat*DX=A'_lat *X+ B'_lat *U -----> DX = A_lat*X +B_lat*U
E_lat         =  zeros(4);
E_lat(1,1)    = -WB.m;
E_lat(2,2)    = -WB.Ixx_b;
E_lat(2,3)    =  WB.Ixz_b;
E_lat(3,2)    =  WB.Ixz_b;
E_lat(3,3)    = -WB.Izz_b;
E_lat(4,4)    =  1;

A_p_lat       =  zeros(4);
A_p_lat(1,1)  =  SD.lat.Yv;
A_p_lat(1,2)  =  SD.lat.Yp;
A_p_lat(1,3)  =  SD.lat.Yr- WB.m*FC.us;
A_p_lat(1,4)  =  WB.m*FC.gravity*cos(RC.Thetas);
A_p_lat(2,1)  =  SD.lat.Lv;
A_p_lat(2,2)  =  SD.lat.Lp;
A_p_lat(2,3)  =  SD.lat.Lr;
A_p_lat(3,1)  =  SD.lat.Nv;
A_p_lat(3,2)  =  SD.lat.Np;
A_p_lat(3,3)  =  SD.lat.Nr;
A_p_lat(4,2)  =  -1;
A_p_lat(4,3)  =  -tan(RC.Thetas);
A_p_lat=-A_p_lat;

B_p_lat       =  zeros(4,2);
B_p_lat(1,2)  = -SD.lat.Ydeltar;
B_p_lat(2,1)  = -SD.lat.Ldeltaa;
B_p_lat(2,2)  = -SD.lat.Ldeltar;
B_p_lat(3,1)  = -SD.lat.Ndeltaa;
B_p_lat(3,2)  = -SD.lat.Ndeltar;

A_lat = E_lat\A_p_lat; %Inv(E_lat)*A_p_lat
B_lat = E_lat\B_p_lat; %Inv(E_lat)*B_p_lat
C = eye(4);
D = zeros(4,2); 

% Funciones de transferencia lateral direccionales;
[TF.lat.Nums_a,TF.lat.Den_a] = ss2tf(A_lat,B_lat,C,D,1);
[TF.lat.Nums_r,TF.lat.Den_r] = ss2tf(A_lat,B_lat,C,D,2);


%% Calculadas las funciones de transferencia guardamos cada una de ellas para trabajar con ellas
% LONGITUDINALES

% Resolvemos la cuartica de estabilidad
eigenvalues_long            =  roots(TF.long.Den);
Cuartica.long.autoval.F     =  eigenvalues_long(3);
Cuartica.long.autoval.SP1    =  eigenvalues_long(1);
Cuartica.long.autoval.SP2    =  eigenvalues_long(2);
Cuartica.long.Omega_n_F     =  abs(Cuartica.long.autoval.F);
Cuartica.long.amort_F       = -(real(Cuartica.long.autoval.F))/Cuartica.long.Omega_n_F ;
Cuartica.long.tau1          = -1/Cuartica.long.autoval.SP1;
Cuartica.long.tau2          = -1/Cuartica.long.autoval.SP2 ;

% Elevador a velocidad
Elevador_veloc.Num       =  TF.long.Nums(1,:);
Elevador_veloc.K         =  polyval(Elevador_veloc.Num,0)/polyval(TF.long.Den,0);
Eigenv_Ele_vel           =  roots(Elevador_veloc.Num);
Elevador_veloc.tau1      =  -1/(Eigenv_Ele_vel(1));
Elevador_veloc.tau2      =  -1/(Eigenv_Ele_vel(2));

% Elevador a �ngulo de Ataque
Elevador_ataque.Num       =  TF.long.Nums(2,:);
Elevador_ataque.K         =  polyval(Elevador_ataque.Num,0)/polyval(TF.long.Den,0);
Eigenv_Ele_ataque         =  roots(Elevador_ataque.Num);
Elevador_ataque.Omega_n   =  abs(Eigenv_Ele_ataque(2));
Elevador_ataque.amort     = -(real(Eigenv_Ele_ataque(2)))/Elevador_ataque.Omega_n;
Elevador_ataque.tau       = -1/Eigenv_Ele_ataque(1);    

% Elevador a �ngulo de Asiento
Elevador_asiento.Num       =  TF.long.Nums(3,:);
Elevador_asiento.K         =  polyval(Elevador_asiento.Num,0)/polyval(TF.long.Den,0);
Eigenv_Ele_asiento         =  roots(Elevador_asiento.Num);
Elevador_asiento.tau1      =  -1/Eigenv_Ele_asiento(1);
Elevador_asiento.tau2      =  -1/Eigenv_Ele_asiento(2);

% Elevador a Velocidad de Cabeceo
Elevador_cabeceo.Num       =  TF.long.Nums(4,:);
Elevador_cabeceo.K         =  polyval(Elevador_cabeceo.Num,0)/polyval(TF.long.Den,0);
Eigenv_Ele_cabeceo        =  roots(Elevador_cabeceo.Num);


%% LATERAL-DIRECCIONALES

% Resolvemos la cuartica de estabilidad
eigenvalues_lateral          = roots(TF.lat.Den_a);
Cuartica.lat.autoval.DR1     = eigenvalues_lateral(1);
Cuartica.lat.autoval.DR2     = eigenvalues_lateral(2);
Cuartica.lat.autoval.CB      = eigenvalues_lateral(3);
Cuartica.lat.autoval.S       = eigenvalues_lateral(4);
Cuartica.lat.Omega_n_DR      = abs(Cuartica.lat.autoval.DR1);
Cuartica.lat.amort_DR        = -(real(Cuartica.lat.autoval.DR1))/Cuartica.lat.Omega_n_DR;
Cuartica.lat.tau_CB          = -1/Cuartica.lat.autoval.CB;
Cuartica.lat.tau_S           = -1/Cuartica.lat.autoval.S;


% Alerones a Velocidad Angular de Balance
Aleron_balance.Num       =  TF.lat.Nums_a(4,:);
Aleron_balance.K         =  polyval(Aleron_balance.Num,0)/polyval(TF.lat.Den_a,0);
Eigenv_balance           =  roots(Aleron_balance.Num);
Aleron_balance.Omega_n   =  abs(Eigenv_balance(1));
Aleron_balance.amort     = -(real(Eigenv_balance(1)))/Aleron_balance.Omega_n; 

% Alerones a Velocidad Angular de Gui�ada
Aleron_vel_guinada.Num      =  TF.lat.Nums_a(3,:);
Aleron_vel_guinada.K        =  polyval(Aleron_vel_guinada.Num,0)/polyval(TF.lat.Den_a,0);
Eigenv_Aleron_vel_guinada   =  roots(Aleron_vel_guinada.Num);
Aleron_vel_guinada.Omega_n  =  abs(Eigenv_Aleron_vel_guinada(1));
Aleron_vel_guinada.amort    = -(real(Eigenv_Aleron_vel_guinada(1)))/Aleron_vel_guinada.Omega_n; 
Aleron_vel_guinada.tau      = -1/Eigenv_Aleron_vel_guinada(3);

% Tim�n de Direcci�n a �ngulo de Resbalamiento
Direccion_resba.Num         =  TF.lat.Nums_r(1,:);
Direccion_resba.K           =  polyval(Direccion_resba.Num,0)/polyval(TF.lat.Den_r,0);
Eigenv_Direccion_resba      =  roots(Direccion_resba.Num );
Direccion_resba.tau1        = -1/Eigenv_Direccion_resba(1);
Direccion_resba.tau2        = -1/Eigenv_Direccion_resba(2);
Direccion_resba.tau3        = -1/Eigenv_Direccion_resba(3);

% Tim�n de Direcci�n a Velocidad Angular de Gui�ada
Direccion_vel_guinada.Num      =  TF.lat.Nums_r(3,:);
Direccion_vel_guinada.K        =  polyval(Direccion_vel_guinada.Num,0)/polyval(TF.lat.Den_r,0);
Eigenv_Direccion_vel_guinada   =  roots(Direccion_vel_guinada.Num);
Direccion_vel_guinada.Omega_n  =  abs(Eigenv_Direccion_vel_guinada(2));
Direccion_vel_guinada.amort    = -(real(Eigenv_Direccion_vel_guinada(2)))/Direccion_vel_guinada.Omega_n; 
Direccion_vel_guinada.tau      = -1/Eigenv_Direccion_vel_guinada(1);

%% FUNCIONES DE TRANSFERENCIA

% Elevador a velocidad
Elev_Vel_K     =  roskam.Elev_vel.K/FC.us;
Elev_Vel_Num1  =  tf([roskam.Elev_vel.tau1 1],1);
Elev_Vel_Num2  =  tf([roskam.Elev_vel.tau2 1],1);
Long_Den1      =  tf([roskam.cuart.tau1_SP 1],1);
Long_Den2      =  tf([roskam.cuart.tau2_SP 1],1);
Long_Den3      =  tf([1/roskam.cuart.w_nF^2 2*roskam.cuart.tsiF/roskam.cuart.w_nF 1],1);
% TF.Elev_vel    =  Elev_Vel_K*Elev_Vel_Num1*Elev_Vel_Num2/(Long_Den1*Long_Den2*Long_Den3);
TF.Elev_vel   =  tf(TF.long.Nums(1,:), TF.long.Den)/FC.us;


% Elevador a ataque
Elev_Ataque_K     =  roskam.Elev_ataque.K;
Elev_Ataque_Num1  =  tf([1/roskam.Elev_ataque.Wn^2 2*roskam.Elev_ataque.Tsi/roskam.Elev_ataque.Wn 1],1);
Elev_Ataque_Num2  =  tf([roskam.Elev_ataque.Tau 1],1);
% TF.Elev_ataque    =  Elev_Ataque_K*Elev_Ataque_Num1*Elev_Ataque_Num2/(Long_Den1*Long_Den2*Long_Den3); 
TF.Elev_ataque    =  tf(TF.long.Nums(2,:), TF.long.Den)/FC.us;

%Elevador a asiento
Elev_Asiento_K     =  roskam.Elev_asiento.K;
Elev_Asiento_Num1  =  tf([roskam.Elev_asiento.Tau1 1],1);
Elev_Asiento_Num2  =  tf([roskam.Elev_asiento.Tau2 1],1);
% TF.Elev_asiento    =  Elev_Asiento_K*Elev_Asiento_Num1*Elev_Asiento_Num2/(Long_Den1*Long_Den2*Long_Den3);
TF.Elev_asiento    =  tf(TF.long.Nums(3,:), TF.long.Den);

%Elevador a cabeceo
Elev_Cabeceo_Num1 =   tf([1 0],1);
% TF.Elev_cabeceo    =  Elev_Asiento_K*Elev_Cabeceo_Num1*Elev_Asiento_Num1*Elev_Asiento_Num2/(Long_Den1*Long_Den2*Long_Den3);
TF.Elev_cabeceo    =  tf(TF.long.Nums(4,:), TF.long.Den);

% Alerones a Velocidad Angular de Balance
Aleron_Balance_K     =  Aleron_balance.K ;
Aleron_Balance_Num1  =  tf([1/Aleron_balance.Omega_n^2 2*Aleron_balance.amort/Aleron_balance.Omega_n 1],1);
Aleron_Balance_Num2  =  tf([1 0],1);
Lat_Den1             =  tf([1/Cuartica.lat.Omega_n_DR^2 2*Cuartica.lat.amort_DR/Cuartica.lat.Omega_n_DR 1],1);
Lat_Den2             =  tf([Cuartica.lat.tau_CB 1],1);
Lat_Den3             =  tf([Cuartica.lat.tau_S 1],1);
TF.Aleron_balance    =  Aleron_Balance_K*Aleron_Balance_Num1*Aleron_Balance_Num2/(Lat_Den1*Lat_Den2*Lat_Den3);

% Alerones a Velocidad Angular de Gui�ada
Aleron_Vel_Guinada_K     =  Aleron_vel_guinada.K ;
Aleron_Vel_Guinada_Num1  =  tf([1/Aleron_vel_guinada.Omega_n^2 2*Aleron_vel_guinada.amort/Aleron_vel_guinada.Omega_n 1],1);
Aleron_Vel_Guinada_Num2  =  tf([Aleron_vel_guinada.tau 1],1);
TF.Aleron_vel_guinada    =  Aleron_Vel_Guinada_K*Aleron_Vel_Guinada_Num1*Aleron_Vel_Guinada_Num2/(Lat_Den1*Lat_Den2*Lat_Den3);

% Tim�n de Direcci�n a �ngulo de Resbalamiento
Direccion_Resbalamiento_K     =  Direccion_resba.K/FC.us ;
Direccion_Resbalamiento_Num1  =  tf([Direccion_resba.tau1 1],1);
Direccion_Resbalamiento_Num2  =  tf([Direccion_resba.tau2 1],1);
Direccion_Resbalamiento_Num3  =  tf([Direccion_resba.tau3 1],1);
TF.Direccion_resbalamiento    =  Direccion_Resbalamiento_K*Direccion_Resbalamiento_Num1*Direccion_Resbalamiento_Num2*Direccion_Resbalamiento_Num3/(Lat_Den1*Lat_Den2*Lat_Den3);


% Tim�n de Direcci�n a Velocidad Angular de Gui�ada
Direccion_Vel_guinada_K       =  Direccion_vel_guinada.K ;
Direccion_Vel_guinada__Num1   =  tf([Direccion_vel_guinada.tau 1],1);
Direccion_Vel_guinada__Num2   =  tf([1/Direccion_vel_guinada.Omega_n^2 2*Direccion_vel_guinada.amort/Direccion_vel_guinada.Omega_n 1],1);
TF.Direccion_Vel_guinada    =  Direccion_Vel_guinada_K*Direccion_Vel_guinada__Num1*Direccion_Vel_guinada__Num2/(Lat_Den1*Lat_Den2*Lat_Den3);


save('utilities/TF_OL.mat','TF')
return






%% DIAGRAMAS DE BODE LONGITUDINALES
%**************************************************************************
%Representaci�n Gr�fica Manual
%**************************************************************************
% vector F.transferencia
TFlongvector = [TF.Elev_vel TF.Elev_ataque TF.Elev_asiento TF.Elev_cabeceo];
% Rango de Frecuencias
OmegaIni = 0.01;
OmegaFin = 100;
Nomega   = 10000;
% Escala logaritmica
vOmega  = logspace(log10(OmegaIni),log10(OmegaFin),Nomega);
% Respuesta en frecuencias
for i=1:4
    titles={'Elevador a Velocidad';'Elevador a \alpha';...
        'Elevador a \theta';'Elevador a q'};
  
    [vMagBodeTemp(i,:) vPhaseBodeTemp(i,:)] = bode(TFlongvector(i),vOmega);
    vMagBodedB(i,:)    =  squeeze(20*log10(vMagBodeTemp(i,:)));
    vPhaseBodeDeg(i,:) =  squeeze(vPhaseBodeTemp(i,:));
    
    
    figure(i);
    set(0,'DefaultLineLineWidth',1.5);
     set(gca,'FontSize',10,'FontName','Times new Roman','box','on')
     grafWidth   = 16;
     grafAR      = 0.7;
     set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     
    subplot(2,1,1);
    plot(vOmega,vMagBodedB(i,:),'Color','b','Linewidth',2)
    grid on; hold all;
    set(gca,'XScale','log');
    set(gca,'Xlim',[OmegaIni OmegaFin]);
     set(0,'DefaultLineLineWidth',1.5);
     set(gca,'FontSize',10,'FontName','Times new Roman','box','on','XScale','log')
     grafWidth   = 16;
     grafAR      = 0.7;
     set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     
    xlabel('Freq (rd/s)');
    ylabel('Amp(dB)');
    title(titles(i));
    subplot(2,1,2);
    semilogx(vOmega,vPhaseBodeDeg(i,:),'Color','r','Linewidth',2)
    grid on; hold all;
    set(gca,'XLim',[OmegaIni OmegaFin]);
     set(0,'DefaultLineLineWidth',1.5);
     set(gca,'FontSize',10,'FontName','Times new Roman','box','on','XScale','log')
     grafWidth   = 16;
     grafAR      = 0.7;
     set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     
    xlabel('Freq (rd/s)');
    ylabel('Phase (Deg)');
    set(gcf,'Color',[1 1 1]);

        format_Grafico = strcat('Diagrama de Bode ',num2str(i));
     format_Grafico = [pwd filesep 'figuras' filesep format_Grafico];

     saveas(gcf,format_Grafico,'epsc');  


end

%% DIAGRAMAS DE BODE LATERALES

%**************************************************************************
%Representaci�n Gr�fica Manual
%**************************************************************************
% vector F.transferencia
TFlatvector = [TF.Aleron_balance TF.Aleron_vel_guinada TF.Direccion_resbalamiento TF.Direccion_Vel_guinada];
% Rango de Frecuencias
OmegaIni = 10^(-5);
OmegaFin = 100;
Nomega   = 10000;
% Escala logaritmica
vOmega  = logspace(log10(OmegaIni),log10(OmegaFin),Nomega);
% Respuesta en frecuencias
for i=1:4
  
    [vMagBodeTemp(i,:) vPhaseBodeTemp(i,:)] = bode(TFlatvector(i),vOmega);
    vMagBodedB(i,:)    =  squeeze(20*log10(vMagBodeTemp(i,:)));
    vPhaseBodeDeg(i,:) =  squeeze(vPhaseBodeTemp(i,:));
    
    titles={'Alerones a p';'Alerones a r';...
       'Tim�n de direcci�n a \beta';'Tim�n de direcci�n a r'};
    
    figure(i+4);
    set(0,'DefaultLineLineWidth',1.5);
    set(gca,'FontSize',10,'FontName','Times new Roman','box','on','XScale','log')
    grafWidth   = 16;
    grafAR      = 0.7;
    set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     
    subplot(2,1,1);
   
    plot(vOmega,vMagBodedB(i,:),'Color','b','Linewidth',2)
    
    grid on; hold all;
    set(gca,'XScale','log');
    set(gca,'Xlim',[OmegaIni OmegaFin]);
    set(0,'DefaultLineLineWidth',1.5);
    set(gca,'FontSize',10,'FontName','Times new Roman','box','on','XScale','log')
    set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     
    xlabel('Freq (rd/s)');
    ylabel('Amp(dB)');
    title(titles(i));
    subplot(2,1,2);
    semilogx(vOmega,vPhaseBodeDeg(i,:),'Color','r','Linewidth',2)
    grid on; hold all;
    set(gca,'XLim',[OmegaIni OmegaFin]);
    set(0,'DefaultLineLineWidth',1.5);
    set(gca,'FontSize',10,'FontName','Times new Roman','box','on','XScale','log')
    set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     
    xlabel('Freq (rd/s)');
    ylabel('Phase (Deg)');
    set(gcf,'Color',[1 1 1]);
    
    format_Grafico = strcat('Diagrama de Bode lat ',num2str(i));
    format_Grafico = [pwd filesep 'figuras' filesep format_Grafico];

    saveas(gcf,format_Grafico,'epsc');  

end

%% RESPUESTA EN FRECUENCIAS
%**************************

% Longitudinal: Obtenci�n de la magnitud a las frecuencias de cada modo propio

vFreqData.long = [roskam.cuart.w_nF, 1/roskam.cuart.tau1_SP, 1/roskam.cuart.tau2_SP];%.*2*pi;
[TF_value]   = freqresp(TFlongvector, vFreqData.long);
TF_value     = squeeze(TF_value);
vMagdB_value.long = 20.*log10(abs(TF_value));

%Lateral-direccional
vFreqData.lat = [Cuartica.lat.Omega_n_DR,  1/Cuartica.lat.tau_CB, -1/Cuartica.lat.tau_S];%.*2*pi;
[TF_value]   = freqresp(TFlatvector, vFreqData.lat);
TF_value     = squeeze(TF_value);
vMagdB_value.lat = 20.*log10(abs(TF_value));

%% PINTAMOS LOS DIAGRAMAS DE NICHOLS LONGITUDINALES

%**************************************************************************
%Representaci�n Gr�fica Manual
%**************************************************************************
% Respuesta en frecuencias
% Rango de Frecuencias
OmegaIni = 0.01;
OmegaFin = 100;
Nomega   = 10000;
% Escala logaritmica
vOmega  = logspace(log10(OmegaIni),log10(OmegaFin),Nomega);
for i=1:4
    titles={'Elevador a Velocidad';'Elevador a \alpha';...
        'Elevador a \theta';'Elevador a q'};
  
    [vMagBodeTemp(i,:) vPhaseBodeTemp(i,:)] = bode(TFlongvector(i),vOmega);
    vMagBodedB(i,:)    =  squeeze(20*log10(vMagBodeTemp(i,:)));
    vPhaseBodeDeg(i,:) =  squeeze(vPhaseBodeTemp(i,:));
    
    figure(i+8);
    plot(vPhaseBodeDeg(i,:),vMagBodedB(i,:),'Color','b','Linewidth',2)
    grid on; hold all;
      set(0,'DefaultLineLineWidth',1.5);
     set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     
    set(gca,'Xlim',[min(vPhaseBodeDeg(i,:)) max(vPhaseBodeDeg(i,:))],'FontSize',10,'FontName','Times new Roman','box','on');
    %xlabel('Freq (rd/s)','FontSize',20);
    xlabel('Phase (Deg)');
    %ylabel('Mag(dB)','FontSize',20)
    ylabel('Mag(dB)')
   % hold on;set(gca,'XTick',[-720:45:720]);
   % hold on;set(gca,'YTick',[-100:10:100]);
    set(gcf,'Color',[1 1 1]);
    title(titles(i));
    
    format_Grafico = strcat('Diagrama de Nichols long ',num2str(i));
    format_Grafico = [pwd filesep 'figuras' filesep format_Grafico];

    saveas(gcf,format_Grafico,'epsc');  
   
    
end

%% PINTAMOS LOS DIAGRAMAS DE NICHOLS LATERALES
%**************************************************************************
%Representaci�n Gr�fica Manual
%**************************************************************************
% Rango de Frecuencias
OmegaIni = 10^(-5);
OmegaFin = 100;
Nomega   = 10000;
% Escala logaritmica
vOmega  = logspace(log10(OmegaIni),log10(OmegaFin),Nomega);
% Respuesta en frecuencias
for i=1:4
    titles={'Alerones a p';'Alerones a r';...
        'Tim�n de direcci�n a \beta';'Tim�n de direcci�n a r'};
  
    [vMagBodeTemp(i,:) vPhaseBodeTemp(i,:)] = bode(TFlatvector(i),vOmega);
    vMagBodedB(i,:)    =  squeeze(20*log10(vMagBodeTemp(i,:)));
    vPhaseBodeDeg(i,:) =  squeeze(vPhaseBodeTemp(i,:));
    
    figure(i+12);
    plot(vPhaseBodeDeg(i,:),vMagBodedB(i,:),'Color','b','Linewidth',2)
    grid on; hold all;
      set(0,'DefaultLineLineWidth',1.5);
     set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     set(gca,'Xlim',[min(vPhaseBodeDeg(i,:)) max(vPhaseBodeDeg(i,:))],'FontSize',10,'FontName','Times new Roman','box','on');
    xlabel('Phase (Deg)');
    ylabel('Mag(dB)')
    hold on;set(gca,'XTick',[-720:45:720]);
    hold on;set(gca,'YTick',[-100:10:100]);
    set(gcf,'Color',[1 1 1]);
    title(titles(i));
    
    format_Grafico = strcat('Diagrama de Nichols lat ',num2str(i));
    format_Grafico = [pwd filesep 'figuras' filesep format_Grafico];
    saveas(gcf,format_Grafico,'epsc');  
end

%% RESPUESTA A ENTRADA ESCALON
%%Longitudinales
%**************************************************************************
%Representaci�n Gr�fica Autom�tica
%**************************************************************************
% for i=1:4
%     
%   figure (i+16)
%   step(TFlongvector(i))
%   set(gcf,'Color',[1 1 1]);
% end

%Manual Estacionario
 TimeEnd=120; %segundos
 StepSize=0.01;
 vTime = [0:StepSize:TimeEnd];
 titles={'Respuesta temporal Velocidad vs Entrada escal�n elevador','Respuesta temporal \alpha vs Entrada escal�n elevador',...
     'Respuesta temporal \theta vs Entrada escal�n elevador', 'Respuesta temporal q vs Entrada escal�n elevador'};
 ylabels={'$$\hat{u}$$','$$\alpha$$ ( $$^{o}$$ )','$$\theta$$ ( $$^{o}$$ )','q (rad/s)'};
 vOut=step(TFlongvector,vTime);
 vEstacionario=[0,0,0,0];
 overshoot_long=[0,0,0,0];
 undershoot_long=[0,0,0,0];
 for i=1:4
    figure(i+16)
    plot(vTime,vOut(:,:,i),'b')
    grid on; hold all;
    set(0,'DefaultLineLineWidth',1.5);
    set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
    set(gca,'Xlim',[0 TimeEnd],'FontSize',10,'FontName','Times new Roman','box','on');
    xlabel('Tiempo(s)');
    ylabel(ylabels(i),'Interpreter','Latex');
    title(titles(i));
    vEstacionario(i)=vOut(length(vTime),:,i);
    plot(vTime,ones(size(vTime))*vEstacionario(i),'--','LineWidth',1);
    
    overshoot_long(i)= -vEstacionario(i)+max(vOut(:,:,i));
    undershoot_long(i)= -vEstacionario(i)+min(vOut(:,:,i));
     
    format_Grafico = strcat('EscalonEstacionario',num2str(i));
    format_Grafico = [pwd filesep 'figuras' filesep format_Grafico];
    saveas(gcf,format_Grafico,'epsc'); 
 end

 
  TimeEnd=13; %segundos
 StepSize=0.01;
 vTime = [0:StepSize:TimeEnd];
 titles={'Transitorio Velocidad vs Entrada escal�n elevador','Transitorio \alpha vs Entrada escal�n elevador',...
     'Transitorio \theta vs Entrada escal�n elevador', 'Transitorio q vs Entrada escal�n elevador'};
 ylabels={'$$\hat{u}$$','$$\alpha$$ ( $$^{o}$$ )','$$\theta$$ ( $$^{o}$$ )','q (rad/s)'};
 vOut=step(TFlongvector,vTime);
%  vEstacionario=[0,0,0,0];
 for i=1:4
    figure(i+20)
    plot(vTime,vOut(:,:,i),'b')
     grid on; hold all;
      set(0,'DefaultLineLineWidth',1.5);
     set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     set(gca,'Xlim',[0 TimeEnd],'FontSize',10,'FontName','Times new Roman','box','on');
    xlabel('Tiempo(s)');
    ylabel(ylabels(i),'Interpreter','Latex');
    title(titles(i));
   % vEstacionario(i)=vOut(length(vTime),:,i);
    % plot(vTime,ones(size(vTime))*vEstacionario(i),'--','LineWidth',1);
        format_Grafico = strcat('EscalonTransitorio',num2str(i));
    format_Grafico = [pwd filesep 'figuras' filesep format_Grafico];
    saveas(gcf,format_Grafico,'epsc'); 
 end

%% RESPUESTA A ENTRADA ESCALON
%%Lateral-Direccional
% 
% for i=1:4
%     
%   figure (i+20)
%   step(TFlatvector(i))
%   set(gcf,'Color',[1 1 1]);
% end

%Manual Estacionario
 TimeEnd=60; %segundos
 StepSize=0.01;
 vTime = [0:StepSize:TimeEnd];
 titles={'Respuesta temporal p vs Entrada escal�n alerones','Respuesta temporal r vs Entrada escal�n alerones',...
     'Respuesta temporal \beta vs Entrada escal�n tim�n de direcci�n', 'Respuesta temporal r vs Entrada escal�n tim�n de direcci�n'};
 ylabels={'p(rad/s)','r(rad/s)','\beta (�)','r(rad/s)'};
 vOut=step(TFlatvector,vTime);
 vEstacionarioLat=[0,0,0,0];
 for i=1:4
    figure(i+24)
    plot(vTime,vOut(:,:,i),'b')
    grid on; hold all;
    set(0,'DefaultLineLineWidth',1.5);
    set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
    set(gca,'Xlim',[0 TimeEnd],'FontSize',10,'FontName','Times new Roman','box','on');
    xlabel('Tiempo(s)');
    ylabel(ylabels(i));
    title(titles(i));
    vEstacionarioLat(i)=vOut(length(vTime),:,i);
    %No ploteo linea horizontal, no existe valor estacionario
    %plot(vTime,ones(size(vTime))*vEstacionarioLat(i),'--','LineWidth',1);
     
    format_Grafico = strcat('EscalonEstacionarioLat',num2str(i));
    format_Grafico = [pwd filesep 'figuras' filesep format_Grafico];
    saveas(gcf,format_Grafico,'epsc'); 
 end
 
  TimeEnd=13; %segundos
 StepSize=0.01;
 vTime = [0:StepSize:TimeEnd];
  titles={'Transitorio p vs Entrada escal�n alerones','Transitorio r vs Entrada escal�n alerones',...
     'Transitorio \beta vs Entrada escal�n tim�n de direcci�n', 'Transitorio r vs Entrada escal�n tim�n de direcci�n'};
 ylabels={'p(rad/s)','r(rad/s)','\beta (�)','r(rad/s)'};
 vOut=step(TFlongvector,vTime);
 for i=1:4
    figure(i+28)
    plot(vTime,vOut(:,:,i),'b')
     grid on; hold all;
      set(0,'DefaultLineLineWidth',1.5);
     set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     set(gca,'Xlim',[0 TimeEnd],'FontSize',10,'FontName','Times new Roman','box','on');
    xlabel('Tiempo(s)');
    ylabel(ylabels(i));
    title(titles(i));
   % vEstacionario(i)=vOut(length(vTime),:,i);
    % plot(vTime,ones(size(vTime))*vEstacionario(i),'--','LineWidth',1);
        format_Grafico = strcat('EscalonTransitorioLat',num2str(i));
    format_Grafico = [pwd filesep 'figuras' filesep format_Grafico];
    saveas(gcf,format_Grafico,'epsc'); 
 end
%% RESPUESTA A ENTRADA ARBITRARIA
%Estacionario

TimeEnd     = 120;
StepSize    = 0.01;
vTime       = [0:StepSize:TimeEnd];
vInpTime    = [0 10 15 TimeEnd];
vInp        = [0 1 1 1];
vInput      = interp1(vInpTime,vInp,vTime);

%Longitudinal estacionaria
for i=1:4
    vTime       = [0:StepSize:TimeEnd];
    vInpTime    = [0 10 15 TimeEnd];
    vInput      = interp1(vInpTime,vInp,vTime);

    vOut_long = lsim(TFlongvector(i),vInput,vTime);

    figure(i+32)
    titles={'Respuesta temporal Velocidad vs rampa elevador','Respuesta temporal \alpha vs rampa elevador',...
         'Respuesta temporal \theta vs rampa elevador', 'Respuesta temporal q vs rampa elevador'};
    ylabels={'u(m/s)','\alpha (�)','\theta (�)','q (rad/s)'};

    plot(vTime,vInput,'r');
    grid on;hold all;
    plot(vTime,vOut_long,'b');
    grid on;hold all;
    title(titles(i))
    set(0,'DefaultLineLineWidth',1.5);
    set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
    set(gca,'Xlim',[0 TimeEnd],'FontSize',10,'FontName','Times new Roman','box','on');
    xlabel('Tiempo(s)');
    ylabel(ylabels(i));
    legend({'Input','Output'})
    set(gcf,'Color',[1 1 1]);
    
    vEstacionario_rampa(i)=vOut_long(length(vTime));
    plot(vTime,ones(size(vTime))*vEstacionario_rampa(i),'--','LineWidth',1);
    
    overshoot_long(i)= -vEstacionario_rampa(i)+max(vOut_long);
    undershoot_long(i)= -vEstacionario_rampa(i)+min(vOut_long);


    format_Grafico = strcat('RampaEstacionario',num2str(i));
    format_Grafico = [pwd filesep 'figuras' filesep format_Grafico];
    saveas(gcf,format_Grafico,'epsc'); 


end
%% Longitudinal transitoria
TimeEnd=20;
for i=1:4
    vTime       = [0:StepSize:TimeEnd];
    vInpTime    = [0 10 15 TimeEnd];
    vInput      = interp1(vInpTime,vInp,vTime);

    vOut_long = lsim(TFlongvector(i),vInput,vTime);

    figure(i+36)
    titles={'Transitorio Velocidad vs rampa elevador','Transitorio \alpha vs rampa elevador',...
         'Transitorio \theta vs rampa elevador', 'Transitorio q vs rampa elevador'};
    ylabels={'u(m/s)','\alpha (�)','\theta (�)','q (rad/s)'};

    plot(vTime,vInput,'r');
    grid on;hold all;
    plot(vTime,vOut_long,'b');
    grid on;hold all;
    title(titles(i))
    set(0,'DefaultLineLineWidth',1.5);
    set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
    set(gca,'Xlim',[0 TimeEnd],'FontSize',10,'FontName','Times new Roman','box','on');
    xlabel('Tiempo(s)');
    ylabel(ylabels(i));
    legend({'Input','Output'})
    set(gcf,'Color',[1 1 1]);
    
    format_Grafico = strcat('RampaTransitorio',num2str(i));
    format_Grafico = [pwd filesep 'figuras' filesep format_Grafico];
    saveas(gcf,format_Grafico,'epsc'); 

end

%% Lateral-direccional estacionaria
TimeEnd=60;
for i=1:4
    vTime       = [0:StepSize:TimeEnd];
    vInpTime    = [0 10 15 TimeEnd];
    vInput      = interp1(vInpTime,vInp,vTime);
    vOut_lat = lsim(TFlatvector(i),vInput,vTime);

    figure(i+40)
    titles={'Respuesta temporal p vs rampa alerones','Respuesta temporal r vs rampa alerones',...
         'Respuesta temporal \beta vs rampa tim�n de direcci�n', 'Respuesta temporal r vs rampa tim�n de direcci�n'};
    ylabels={'p(rad/s)','r(rad/s)','\beta (�)','r(rad/s)'};


    plot(vTime,vInput,'r');
    grid on;hold all;
    plot(vTime,vOut_lat,'b');
    grid on;hold all;
    title(titles(i))
        set(0,'DefaultLineLineWidth',1.5);
        set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
        set(gca,'Xlim',[0 TimeEnd],'FontSize',10,'FontName','Times new Roman','box','on');
        xlabel('Tiempo(s)');
        ylabel(ylabels(i));
    legend({'Input','Output'})
    set(gcf,'Color',[1 1 1]);
    
    format_Grafico = strcat('RampaEstacionarioLat',num2str(i));
    format_Grafico = [pwd filesep 'figuras' filesep format_Grafico];
    saveas(gcf,format_Grafico,'epsc'); 


end
%Lateral-direccional transitoria
TimeEnd=16;
for i=1:4
    vTime       = [0:StepSize:TimeEnd];
    vInpTime    = [0 10 15 TimeEnd];
    vInput      = interp1(vInpTime,vInp,vTime);
    vOut_lat = lsim(TFlatvector(i),vInput,vTime);

    figure(i+44)
    titles={'Transitorio p vs rampa alerones','Transitorio r vs rampa alerones',...
         'Transitorio \beta vs rampa tim�n de direcci�n', 'Transitorio r vs rampa tim�n de direcci�n'};
    ylabels={'p(rad/s)','r(rad/s)','\beta (�)','r(rad/s)'};


    plot(vTime,vInput,'r');
    grid on;hold all;
    plot(vTime,vOut_lat,'b');
    grid on;hold all;
    title(titles(i))
        set(0,'DefaultLineLineWidth',1.5);
        set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
        set(gca,'Xlim',[0 TimeEnd],'FontSize',10,'FontName','Times new Roman','box','on');
        xlabel('Tiempo(s)');
        ylabel(ylabels(i));
    legend({'Input','Output'})
    set(gcf,'Color',[1 1 1]);
    
    format_Grafico = strcat('RampaTransitorioLat',num2str(i));
    format_Grafico = [pwd filesep 'figuras' filesep format_Grafico];
    saveas(gcf,format_Grafico,'epsc'); 

end
