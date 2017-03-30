%%
clc
close all
clear all
%% CONVERSIÓN DE UNIDADES A S.I.
% Definicion de los parametros de la practica:Cessna 310 Approach configuration
% Flight Condition FC
FC.hs  =   0;           % Altura en [pies]
FC.Ms  =   0.124;       % Número de Mach
FC.us  =   137.9;       % Velocidad condicion de referencia [pies/segundos]
FC.qs  =   22.6 ;       % presion dinámica condicion de referencia [libra/pies^2] 
FC.alfa_xbxs = 6.6;     % ángulo entre la xB y xs en [º]
%Cambiamos las unidades
ft2m   =   1/3.28084;      %Factor de conversion pies a metros 
lb2kg  =   0.453592;       %Factor de conversion libra a kg
lbft2Pa =  47.88;          %Factor de conversion presion dinámica a pascales
slgft2kgm = 1.355817962;   %Factor de conversion slug*pie^2 a kg*m^2
deg2rad   = pi/180;        %Factor de conversion grados a radianes

% Flight Condition FC
FC.hs  =   FC.hs*ft2m;                    % Altura en [m]
FC.us  =   FC.us*ft2m;                    % Velocidad condicion de referencia [m/segundos]
FC.qs  =   FC.qs*lbft2Pa ;                % presion dinámica condicion de referencia [libra/pies^2] 
FC.alfa_xbxs = FC.alfa_xbxs*deg2rad ;     % ángulo entre la xB y xs en [º]
 
% Datos Geometricos GEO
Geo.Sw  =  175;          % Superficie alar [pies^2]
Geo.c   =  4.79;         % Cuerda [pies]
Geo.b   =  36.9;         % envergadura [pies]

Geo.Sw  =  Geo.Sw*ft2m^2;      % Superficie alar [metros^2]
Geo.c   =  Geo.c*ft2m;         % Cuerda [metros]
Geo.b   =  Geo.b*ft2m;         % envergadura [metros]

% Pesos y Equilibrado WB(weight and balance)

WB.w     =  4600;        % masa en [libras]
WB.Ixx_b =  8884;        % Ixx en ejes cuerpo [slgft^2]/[lbm·ft2]
WB.Iyy_b =  1939;        % Ixx en ejes cuerpo [slgft^2]/[lbm·ft2]
WB.Izz_b =  11001;       % Ixx en ejes cuerpo [slgft^2]/[lbm·ft2]
WB.Ixz_b =  0;           % Ixx en ejes cuerpo [slgft^2]/[lbm·ft2]

WB.w     = WB.w*lb2kg ;                % masa en [libras]
WB.Ixx_b = WB.Ixx_b*slgft2kgm ;        % Ixx en ejes cuerpo [kg*m^2]
WB.Iyy_b = WB.Iyy_b*slgft2kgm ;        % Ixx en ejes cuerpo [kg*m^2]
WB.Izz_b = WB.Izz_b*slgft2kgm ;        % Ixx en ejes cuerpo [kg*m^2]
WB.Ixz_b = WB.Ixz_b*slgft2kgm ;        % Ixx en ejes cuerpo [kg*m^2]


% Condicion de Referencia RF

RF.Cls   = 1.163;   %Coeficiente de sustentacion referencia
RF.Cds   = 0.1710;  %Coeficiente de resistencia referencia
RF.Ct_xs = 0.1710;  %Coeficiente de traccion referencia en eje Xs
RF.Cms   = 0;       %Coeficiente de momento referencia
RF.Cmts  = 0;       %Coeficiente de momento de traccion referencia

% Derivadas de Estabilidad Adimensionales longitudimales ASD (Stability Derivatives)

ASD.long.CD_0      =   0.0974;
ASD.long.CD_u      =   0;
ASD.long.CD_alfa   =   0.650;
ASD.long.CT_xu     =  -0.513;
ASD.long.CT_alfa   =   0;
ASD.long.CL_0      =   0.640;
ASD.long.CL_u      =   0;
ASD.long.CL_alfa   =   4.58;
ASD.long.CL_alfap  =   4.1;  % CL respecto a alfa punto adimensional
ASD.long.CL_q      =   8.4;
ASD.long.Cm_0      =   0.1;
ASD.long.Cm_u      =   0;
ASD.long.Cm_alfa   =  -0.619;
ASD.long.Cm_alfap  =  -11.40;
ASD.long.Cm_q      =  -25.1;
ASD.long.Cm_Tu     =   0;
ASD.long.Cm_Talfa  =   0;
ASD.long.CD_deltae =   0;
ASD.long.CL_deltae =   0.77;
ASD.long.Cm_deltae =  -2.16;

% Derivadas de Estabilidad Adimensionales lateral direccionales ASD (Stability Derivatives)

ASD.lat.Cl_beta   =  -0.0965;
ASD.lat.Cl_p      =  -0.566;
ASD.lat.Cl_r      =   0.2433;
ASD.lat.CY_beta   =  -0.0965;
ASD.lat.CY_p      =  -0.2897;
ASD.lat.CY_r      =   0.355;
ASD.lat.Cn_beta   =   0.1683;
ASD.lat.Cn_T_beta =   0;
ASD.lat.Cn_p      =  -0.1021;
ASD.lat.Cn_r      =  -0.1947;
ASD.lat.Cl_deltaa =   0.1720;
ASD.lat.Cl_deltar =  -0.0192;
ASD.lat.CY_deltaa =   0;
ASD.lat.CY_deltar =  -0.230;
ASD.lat.Cn_deltaa =  -0.0676;
ASD.lat.Cn_deltar =   0.1152;

%% Una vez se tienen los datos, se pasa a obtener las derivadas de estabilidad dimensionales
% Derivadas de estabilidad dimensionales

[T, a, P, rho] = atmosisa(FC.hs);

%Adimensionalizadores longitudinales

Adimlong_1 = rho*FC.us*Geo.Sw;             %rho*Us*S
Adimlong_2 = rho*FC.us*Geo.Sw/2;           %(rho*Us*S)/2
Adimlong_3 = rho*(FC.us)^2*Geo.Sw/2;       %(rho*(Us)^2*S)/2
Adimlong_4 = rho*Geo.Sw*Geo.c/4;           %(rho*S*c)/4
Adimlong_5 = rho*FC.us*Geo.Sw*Geo.c/4;     %(rho*Us*S*c)/4
Adimlong_6 = rho*FC.us*Geo.Sw*Geo.c/2;     %(rho*Us*S*c)/2
Adimlong_7 = rho*Geo.Sw*(Geo.c)^2/4;       %(rho*S*c^2)/4
Adimlong_8 = rho*FC.us*Geo.Sw*(Geo.c)^2/4; %(rho*Us*S*c^2)/4
Adimlong_9 = rho*(FC.us)^2*Geo.Sw*Geo.c/2; %(rho*(Us)^2*S*c)/2
Adimlong_10= rho*(FC.us)^2*Geo.Sw*Geo.c/4; %(rho*(Us)^2*S*c)/4

%Derivadas de Estabilidad ASD longitudinales a partir de las que se tienen
%Hipotesis de alfa pequeña 

C_xs     =    RF.Ct_xs-RF.Cds;
C_zs     =  - RF.Cls;  %RF.Ct_xs
C_ms     =    RF.Cms - RF.Cmts;

C_xu     =  ASD.long.CT_xu - ASD.long.CD_u;
C_zu     = -ASD.long.CL_u ;
C_mu     =  ASD.long.Cm_u- ASD.long.Cm_Tu ;

C_xalfa  =  ASD.long.CT_alfa+RF.Cls - ASD.long.CD_alfa; %%
C_zalfa  = -ASD.long.CL_alfa; %%
C_malfa  =  ASD.long.Cm_alfa - ASD.long.Cm_Talfa; %%

C_zalfap = -ASD.long.CL_alfap; %%
C_malfap =  ASD.long.Cm_alfap; %%

C_zq     = -ASD.long.CL_q;
C_mq     =  ASD.long.Cm_q;

C_xdeltae = - ASD.long.CD_deltae;
C_zdeltae = -ASD.long.CL_deltae;
C_mdeltae =  ASD.long.Cm_deltae;
  

SD.long.Xu      = Adimlong_1*C_xs + Adimlong_2*C_xu; %X_u 
SD.long.Zu      = Adimlong_1*C_zs + Adimlong_2*C_zu; %Z_u
SD.long.Mu      = Adimlong_6*C_mu;                   %M_u
SD.long.Xw      = Adimlong_2*C_xalfa;                %X_w
SD.long.Zw      = Adimlong_2*C_zalfa;                %Z_w
SD.long.Mw      = Adimlong_6*C_malfa;                %M_w
SD.long.Zwp     = Adimlong_4*C_zalfap;               %Z_wp
SD.long.Mwp     = Adimlong_7*C_malfap;               %M_wp
SD.long.Zq      = Adimlong_5*C_zq;                   %Z_q
SD.long.Mq      = Adimlong_8*C_mq;                   %M_q
SD.long.Xdeltae = Adimlong_3*C_xdeltae;              %X_deltae
SD.long.Zdeltae = Adimlong_3*C_zdeltae;              %Z_deltae
SD.long.Mdeltae = Adimlong_9*C_mdeltae;              %M_deltae

%Adimensionalizadores lateral-direccionales

Adimlat_1 = rho*FC.us*Geo.Sw/2;              %(rho*Us*S)/2
Adimlat_2 = rho*(FC.us)*Geo.Sw*Geo.b/2;       %(rho*Us*S*b)/2
Adimlat_3 = rho*(FC.us)*Geo.Sw*Geo.b/4;       %(rho*Us*S*b)/4
Adimlat_4 = rho*(FC.us)*Geo.Sw*(Geo.b)^2/4;   %(rho*Us*S*b^2)/4
Adimlat_5 = rho*(FC.us)^2*Geo.Sw*Geo.b/2;     %(rho*Us^2*S*b)/2
Adimlat_6 = rho*(FC.us)^2*Geo.Sw/2;           %(rho*Us^2*S)/2

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

%% SISTEMA ROSKAM

roskam.X_u      = -FC.qs*Geo.Sw*(ASD.long.CD_u+2*RF.Cds)/(WB.w*FC.us);
roskam.X_Tu     =  FC.qs*Geo.Sw*(ASD.long.CT_xu+2*RF.Ct_xs)/(WB.w*FC.us);
roskam.X_alfa   = -FC.qs*Geo.Sw*(ASD.long.CD_alfa-RF.Cls)/WB.w;
roskam.X_deltae = -FC.qs*Geo.Sw*ASD.long.CD_deltae/WB.w;
roskam.Z_u      = -FC.qs*Geo.Sw*(ASD.long.CL_u+2*RF.Cls)/(WB.w*FC.us);
roskam.Z_alfa   = -FC.qs*Geo.Sw*(ASD.long.CL_alfa-RF.Cds)/WB.w;
roskam.Z_alfap  = -FC.qs*Geo.Sw*Geo.c*ASD.long.CL_alfap/(2*WB.w*FC.us);
roskam.Z_q      = -FC.qs*Geo.Sw*Geo.c*ASD.long.CL_q/(2*WB.w*FC.us);
roskam.Z_deltae = -FC.qs*Geo.Sw*ASD.long.CL_deltae/WB.w;
roskam.M_u      =  FC.qs*Geo.Sw*Geo.c*(ASD.long.Cm_u+2*RF.Cms)/(WB.Iyy_b*FC.us);
roskam.M_Tu     =  FC.qs*Geo.Sw*Geo.c*(ASD.long.Cm_Tu+2*RF.Cmts)/(WB.Iyy_b*FC.us);
roskam.M_alfa   =  FC.qs*Geo.Sw*Geo.c*ASD.long.Cm_alfa/WB.Iyy_b;
roskam.M_Talfa  =  FC.qs*Geo.Sw*Geo.c*ASD.long.Cm_Talfa/WB.Iyy_b;
roskam.M_alfap  =  FC.qs*Geo.Sw*Geo.c^2*ASD.long.Cm_alfap/(2*WB.Iyy_b*FC.us);
roskam.M_q      =  FC.qs*Geo.Sw*Geo.c^2*ASD.long.Cm_q/(2*WB.Iyy_b*FC.us);
roskam.M_deltae =  FC.qs*Geo.Sw*Geo.c*ASD.long.Cm_deltae/WB.Iyy_b;
g               = 9.81;
thetas          =  0; 
syms s

%Matriz MS

MSS = zeros (3);
MS = zeros (3);
ML = zeros (3);

MSS(3,3) = 1;

MS(1,1) = 1;
MS(2,2) = (FC.us-roskam.Z_alfap);
MS(2,3) = -(roskam.Z_q+FC.us);
MS(3,2) = -roskam.M_alfap;
MS(3,3) = -roskam.M_q;

ML(1,1) = -roskam.X_u-roskam.X_Tu;
ML(1,2) = -roskam.X_alfa;
ML(1,3) =  g*cos(thetas);
ML(2,1) = -roskam.Z_u;
ML(2,2) = -roskam.Z_alfa;
ML(2,3) = g*sin(thetas);
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


%%CRAMER de cada modo

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


%% Metodo adimensional Longitudinal
%Construimos matriz con las ecuaciones adimensionales dimensionalizada y hecha su transformada de Laplace
syms s

%Matriz MS

MSS = zeros (4);
MS = zeros (4);
ML = zeros (4);
MS(1,1) =  (2*WB.w/(rho*Geo.Sw*FC.us^2));
ML(1,1) = -(C_xu/FC.us);
ML(1,2) = -C_xalfa;%%
ML(1,3) = -C_zs;
ML(2,1) = -((C_zu+2*C_zs)/FC.us);
MS(2,2) =  (((4*WB.w/(rho*Geo.Sw*Geo.c))-C_zalfap)*Geo.c)/(2*FC.us);%%
ML(2,2) = -C_zalfa;%%
MS(2,3) = -(((4*WB.w/(rho*Geo.Sw*Geo.c))+C_zq)*Geo.c)/(2*FC.us);
ML(3,1) = -C_mu/FC.us;
ML(3,2) = -C_malfa;%%
MS(3,2) = -C_malfap*Geo.c/(2*FC.us);%%
MS(3,3) = -C_mq*Geo.c/(2*FC.us);
MS(3,4) =  WB.Iyy_b/(rho*Geo.Sw*(Geo.c/2))*(1/FC.us^2);
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

%%CRAMER de cada modo

% Elevador a velocidad
K  = MS;
KK = ML;
K(:,1) =0 ;
KK(1,1)= C_xdeltae;
KK(2,1)= C_zdeltae;
KK(3,1)= C_mdeltae;


Cramer.Elev_vel.sist    = simplify(s^2*MSS+s*K+KK);
Cramer.Elev_vel.ec      = det(Cramer.Elev_vel.sist);
Cramer.Elev_vel.coef    = sym2poly(Cramer.Elev_vel.ec);
Cramer.Elev_vel.raices  = roots(Cramer.Elev_vel.coef);


% Elevador a Ángulo de Ataque

K  = MS;
KK = ML;
K(:,2) =0 ;
KK(1,2)= C_xdeltae;
KK(2,2)= C_zdeltae;
KK(3,2)= C_mdeltae;


Cramer.Elev_ataque.sist    = simplify(s^2*MSS+s*K+KK);
Cramer.Elev_ataque.ec      = det(Cramer.Elev_ataque.sist);
Cramer.Elev_ataque.coef    = sym2poly(Cramer.Elev_ataque.ec);
Cramer.Elev_ataque.raices  = roots(Cramer.Elev_ataque.coef);
Cramer.Elev_ataque.Wn      = abs(Cramer.Elev_ataque.raices(2));
Cramer.Elev_ataque.Tsi     =-real(Cramer.Elev_ataque.raices(2))/Cramer.Elev_ataque.Wn;
Cramer.Elev_ataque.Tau     = -1/Cramer.Elev_ataque.raices(1);
Cramer.Elev_ataque.K       = -Cramer.Elev_ataque.coef(1)*Cramer.Elev_ataque.Wn^2*(-1/Cramer.Elev_ataque.Tau)/((-1/Cramer.cuart.tau1_SP)*(-1/Cramer.cuart.tau2_SP)*Cramer.cuart.w_nF^2*Cramer.cuart.Cuarticalong_coef(1));
% Elevador a Ángulo de Asiento
K  = MS;
KK = ML;
KKK =MSS;
K(:,3) =0 ;
KK(1,3)= C_xdeltae;
KK(2,3)= C_zdeltae;
KK(3,3)= C_mdeltae;
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
E(1,1)    = -WB.w;
E(2,2)    =  SD.long.Zwp - WB.w;
E(2,3)    =  0;
E(3,2)    =  SD.long.Mwp;
E(3,3)    =  SD.long.Mq;
E(3,4)    = -WB.Iyy_b;
E(4,3)    =  1;

Theta_s   =  0;        %theta_s theta referencia 
g         =  9.81;      %gravedad m/s^2
A_p       =  zeros(4);
A_p(1,1)  =  SD.long.Xu;
A_p(1,2)  =  SD.long.Xw;
A_p(1,3)  = -WB.w*g*cos(Theta_s);
A_p(2,1)  =  SD.long.Zu;
A_p(2,2)  =  SD.long.Zw;
A_p(2,3)  = -WB.w*g*sin(Theta_s);
A_p(2,4)  =  SD.long.Zq + WB.w*FC.us;
A_p(3,1)  =  SD.long.Mu;
A_p(3,2)  =  SD.long.Mw;
A_p(3,4)  =  0;
A_p(4,4)  =  -1;
A_p=-A_p;

B_p       =  zeros(4,1);
B_p(1,1)  = -SD.long.Xdeltae;
B_p(2,1)  = -SD.long.Zdeltae;
B_p(3,1)  = -SD.long.Mdeltae;

E_inv = inv(E);

A = E_inv*A_p;
B = E_inv*B_p;
C = eye(4);
D = zeros(4,1); 
% Funciones de transferencia longitudinales


[TF.long.Nums,TF.long.Den] = ss2tf(A,B,C,D);


%% Ahora el lateral-direccional
%matriz E_lat*DX=A'_lat *X+ B'_lat *U -----> DX = A_lat*X +B_lat*U
E_lat         =  zeros(4);
E_lat(1,1)    = -WB.w;
E_lat(2,2)    = -WB.Ixx_b;
E_lat(2,3)    =  WB.Ixz_b;
E_lat(3,2)    =  WB.Ixz_b;
E_lat(3,3)    = -WB.Izz_b;
E_lat(4,4)    =  1;

A_p_lat       =  zeros(4);
A_p_lat(1,1)  =  SD.lat.Yv;
A_p_lat(1,2)  =  SD.lat.Yp;
A_p_lat(1,3)  =  SD.lat.Yr- WB.w*FC.us;
A_p_lat(1,4)  =  WB.w*g*cos(Theta_s);
A_p_lat(2,1)  =  SD.lat.Lv;
A_p_lat(2,2)  =  SD.lat.Lp;
A_p_lat(2,3)  =  SD.lat.Lr;
A_p_lat(3,1)  =  SD.lat.Nv;
A_p_lat(3,2)  =  SD.lat.Np;
A_p_lat(3,3)  =  SD.lat.Nr;
A_p_lat(4,2)  =  -1;
A_p_lat(4,3)  =  -tan(Theta_s);
A_p_lat=-A_p_lat;

B_p_lat       =  zeros(4,2);
B_p_lat(1,2)  = -SD.lat.Ydeltar;
B_p_lat(2,1)  = -SD.lat.Ldeltaa;
B_p_lat(2,2)  = -SD.lat.Ldeltar;
B_p_lat(3,1)  = -SD.lat.Ndeltaa;
B_p_lat(3,2)  = -SD.lat.Ndeltar;

E_inv_lat = inv(E_lat);

A_lat = E_inv_lat*A_p_lat;
B_lat = E_inv_lat*B_p_lat;
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
Cuartica.long.tau1          = -1/ Cuartica.long.autoval.SP1;
Cuartica.long.tau2          = -1/Cuartica.long.autoval.SP2 ;

% Elevador a velocidad

Elevador_veloc.Num       =  TF.long.Nums(1,:);
Elevador_veloc.K         =  polyval(Elevador_veloc.Num,0)/polyval(TF.long.Den,0);
Eigenv_Ele_vel           =  roots(Elevador_veloc.Num);
Elevador_veloc.tau1      =  -1/(Eigenv_Ele_vel(1));
Elevador_veloc.tau2      =  -1/(Eigenv_Ele_vel(2));


% Elevador a Ángulo de Ataque

Elevador_ataque.Num       =  TF.long.Nums(2,:);
Elevador_ataque.K         =  polyval(Elevador_ataque.Num,0)/polyval(TF.long.Den,0);
Eigenv_Ele_ataque         =  roots(Elevador_ataque.Num);
Elevador_ataque.Omega_n   =  abs(Eigenv_Ele_ataque(2));
Elevador_ataque.amort     = -(real(Eigenv_Ele_ataque(2)))/Elevador_ataque.Omega_n;
Elevador_ataque.tau       = -1/Eigenv_Ele_ataque(1);    

% Elevador a Ángulo de Asiento

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

% Alerones a Velocidad Angular de Guiñada

Aleron_vel_guinada.Num      =  TF.lat.Nums_a(3,:);
Aleron_vel_guinada.K        =  polyval(Aleron_vel_guinada.Num,0)/polyval(TF.lat.Den_a,0);
Eigenv_Aleron_vel_guinada   =  roots(Aleron_vel_guinada.Num);
Aleron_vel_guinada.Omega_n  =  abs(Eigenv_Aleron_vel_guinada(1));
Aleron_vel_guinada.amort    = -(real(Eigenv_Aleron_vel_guinada(1)))/Aleron_vel_guinada.Omega_n; 
Aleron_vel_guinada.tau      = -1/Eigenv_Aleron_vel_guinada(3);

% Timón de Dirección a Ángulo de Resbalamiento

Direccion_resba.Num         =  TF.lat.Nums_r(1,:);
Direccion_resba.K           =  polyval(Direccion_resba.Num,0)/polyval(TF.lat.Den_r,0);
Eigenv_Direccion_resba      =  roots(Direccion_resba.Num );
Direccion_resba.tau1        = -1/Eigenv_Direccion_resba(1);
Direccion_resba.tau2        = -1/Eigenv_Direccion_resba(2);
Direccion_resba.tau3        = -1/Eigenv_Direccion_resba(3);

% Timón de Dirección a Velocidad Angular de Guiñada

Direccion_vel_guinada.Num      =  TF.lat.Nums_r(3,:);
Direccion_vel_guinada.K        =  polyval(Direccion_vel_guinada.Num,0)/polyval(TF.lat.Den_r,0);
Eigenv_Direccion_vel_guinada   =  roots(Direccion_vel_guinada.Num);
Direccion_vel_guinada.Omega_n  =  abs(Eigenv_Direccion_vel_guinada(2));
Direccion_vel_guinada.amort    = -(real(Eigenv_Direccion_vel_guinada(2)))/Direccion_vel_guinada.Omega_n; 
Direccion_vel_guinada.tau      = -1/Eigenv_Direccion_vel_guinada(1);

%% A partir de los analisis realizados construimos las funciones transferencia

% Elevador a velocidad

Elev_Vel_K     =  roskam.Elev_vel.K/FC.us;
Elev_Vel_Num1  =  tf([roskam.Elev_vel.tau1 1],[1]);
Elev_Vel_Num2  =  tf([roskam.Elev_vel.tau2 1],[1]);
Long_Den1      =  tf([roskam.cuart.tau1_SP 1],[1]);
Long_Den2      =  tf([roskam.cuart.tau2_SP 1],[1]);
Long_Den3      =  tf([1/roskam.cuart.w_nF^2 2*roskam.cuart.tsiF/roskam.cuart.w_nF 1],[1]);
TF.Elev_vel    =  Elev_Vel_K*Elev_Vel_Num1*Elev_Vel_Num2/(Long_Den1*Long_Den2*Long_Den3);

% Elevador a ataque

Elev_Ataque_K     =  roskam.Elev_ataque.K;
Elev_Ataque_Num1  =  tf([1/roskam.Elev_ataque.Wn^2 2*roskam.Elev_ataque.Tsi/roskam.Elev_ataque.Wn 1],[1]);
Elev_Ataque_Num2  =  tf([roskam.Elev_ataque.Tau 1],[1]);
TF.Elev_ataque    =  Elev_Ataque_K*Elev_Ataque_Num1*Elev_Ataque_Num2/(Long_Den1*Long_Den2*Long_Den3); 

%Elevador a asiento

Elev_Asiento_K     =  roskam.Elev_asiento.K;
Elev_Asiento_Num1  =  tf([roskam.Elev_asiento.Tau1 1],[1]);
Elev_Asiento_Num2  =  tf([roskam.Elev_asiento.Tau2 1],[1]);
TF.Elev_asiento    =  Elev_Asiento_K*Elev_Asiento_Num1*Elev_Asiento_Num2/(Long_Den1*Long_Den2*Long_Den3);

%Elevador a cabeceo

Elev_Cabeceo_Num1 =   tf([1 0],[1]);
TF.Elev_cabeceo    =  Elev_Asiento_K*Elev_Cabeceo_Num1*Elev_Asiento_Num1*Elev_Asiento_Num2/(Long_Den1*Long_Den2*Long_Den3);


% Alerones a Velocidad Angular de Balance

Aleron_Balance_K     =  Aleron_balance.K ;
Aleron_Balance_Num1  =  tf([1/Aleron_balance.Omega_n^2 2*Aleron_balance.amort/Aleron_balance.Omega_n 1],[1]);
Aleron_Balance_Num2  =  tf([1 0],[1]);
Lat_Den1             =  tf([1/Cuartica.lat.Omega_n_DR^2 2*Cuartica.lat.amort_DR/Cuartica.lat.Omega_n_DR 1],[1]);
Lat_Den2             =  tf([Cuartica.lat.tau_CB 1],[1]);
Lat_Den3             =  tf([Cuartica.lat.tau_S 1],[1]);
TF.Aleron_balance    =  Aleron_Balance_K*Aleron_Balance_Num1*Aleron_Balance_Num2/(Lat_Den1*Lat_Den2*Lat_Den3);

% Alerones a Velocidad Angular de Guiñada

Aleron_Vel_Guinada_K     =  Aleron_vel_guinada.K ;
Aleron_Vel_Guinada_Num1  =  tf([1/Aleron_vel_guinada.Omega_n^2 2*Aleron_vel_guinada.amort/Aleron_vel_guinada.Omega_n 1],[1]);
Aleron_Vel_Guinada_Num2  =  tf([Aleron_vel_guinada.tau 1],[1]);
TF.Aleron_vel_guinada    =  Aleron_Vel_Guinada_K*Aleron_Vel_Guinada_Num1*Aleron_Vel_Guinada_Num2/(Lat_Den1*Lat_Den2*Lat_Den3);

% Timón de Dirección a Ángulo de Resbalamiento

Direccion_Resbalamiento_K     =  Direccion_resba.K/FC.us ;
Direccion_Resbalamiento_Num1  =  tf([Direccion_resba.tau1 1],[1]);
Direccion_Resbalamiento_Num2  =  tf([Direccion_resba.tau2 1],[1]);
Direccion_Resbalamiento_Num3  =  tf([Direccion_resba.tau3 1],[1]);
TF.Direccion_resbalamiento    =  Direccion_Resbalamiento_K*Direccion_Resbalamiento_Num1*Direccion_Resbalamiento_Num2*Direccion_Resbalamiento_Num3/(Lat_Den1*Lat_Den2*Lat_Den3);


% Timón de Dirección a Velocidad Angular de Guiñada

Direccion_Vel_guinada_K       =  Direccion_vel_guinada.K ;
Direccion_Vel_guinada__Num1   =  tf([Direccion_vel_guinada.tau 1],[1]);
Direccion_Vel_guinada__Num2   =  tf([1/Direccion_vel_guinada.Omega_n^2 2*Direccion_vel_guinada.amort/Direccion_vel_guinada.Omega_n 1],[1]);


TF.Direccion_Vel_guinada    =  Direccion_Vel_guinada_K*Direccion_Vel_guinada__Num1*Direccion_Vel_guinada__Num2/(Lat_Den1*Lat_Den2*Lat_Den3);




%% DIAGRAMAS DE BODE LONGITUDINALES

%**************************************************************************
%Representación Gráfica Manual
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
%Representación Gráfica Manual
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
       'Timón de dirección a \beta';'Timón de dirección a r'};
    
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

% Longitudinal: Obtención de la magnitud a las frecuencias de cada modo propio

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
%Representación Gráfica Manual
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
%Representación Gráfica Manual
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
        'Timón de dirección a \beta';'Timón de dirección a r'};
  
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
%Representación Gráfica Automática
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
 titles={'Respuesta temporal Velocidad vs Entrada escalón elevador','Respuesta temporal \alpha vs Entrada escalón elevador',...
     'Respuesta temporal \theta vs Entrada escalón elevador', 'Respuesta temporal q vs Entrada escalón elevador'};
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
 titles={'Transitorio Velocidad vs Entrada escalón elevador','Transitorio \alpha vs Entrada escalón elevador',...
     'Transitorio \theta vs Entrada escalón elevador', 'Transitorio q vs Entrada escalón elevador'};
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
 titles={'Respuesta temporal p vs Entrada escalón alerones','Respuesta temporal r vs Entrada escalón alerones',...
     'Respuesta temporal \beta vs Entrada escalón timón de dirección', 'Respuesta temporal r vs Entrada escalón timón de dirección'};
 ylabels={'p(rad/s)','r(rad/s)','\beta (º)','r(rad/s)'};
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
  titles={'Transitorio p vs Entrada escalón alerones','Transitorio r vs Entrada escalón alerones',...
     'Transitorio \beta vs Entrada escalón timón de dirección', 'Transitorio r vs Entrada escalón timón de dirección'};
 ylabels={'p(rad/s)','r(rad/s)','\beta (º)','r(rad/s)'};
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
    ylabels={'u(m/s)','\alpha (º)','\theta (º)','q (rad/s)'};

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
    ylabels={'u(m/s)','\alpha (º)','\theta (º)','q (rad/s)'};

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
         'Respuesta temporal \beta vs rampa timón de dirección', 'Respuesta temporal r vs rampa timón de dirección'};
    ylabels={'p(rad/s)','r(rad/s)','\beta (º)','r(rad/s)'};


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
         'Transitorio \beta vs rampa timón de dirección', 'Transitorio r vs rampa timón de dirección'};
    ylabels={'p(rad/s)','r(rad/s)','\beta (º)','r(rad/s)'};


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
