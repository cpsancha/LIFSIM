
clc
close all
clear all
%% HAY QUE CAMBIAR UNIDADES
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

%% Metodo adimensional Longitudinal
%Construimos matriz con las ecuaciones adimensionales dimensionalizada y hecha su transformada de Laplace
syms s

%Matriz MS

MSS = zeros (2);
MS = zeros (2);
ML = zeros (2);

C_zalfap=0;
C_zq=0;

MS(1,1) =  (((4*WB.w/(rho*Geo.Sw*Geo.c))-C_zalfap)*Geo.c)/(2*FC.us);%%
ML(1,1) = -C_zalfa;%%
ML(1,2) = -(((4*WB.w/(rho*Geo.Sw*Geo.c))+C_zq)*Geo.c)/(2*FC.us);
ML(2,1) = C_malfa;%%
MS(2,1) = C_malfap*Geo.c/(2*FC.us);%%
ML(2,2) = C_mq*Geo.c/(2*FC.us);%%%%% Cmq
MS(2,2) =  -WB.Iyy_b/(rho*Geo.Sw*(Geo.c/2))*(1/FC.us^2);


Cramer.cuart.matriz            = s^2*MSS+s*MS+ML;
Cramer.cuart.deter_matriz      = det(Cramer.cuart.matriz);
Cramer.cuart.Cuarticalong_coef = sym2poly(Cramer.cuart.deter_matriz); 
Cramer.cuart.Autovalores       = roots(Cramer.cuart.Cuarticalong_coef);
Cramer.cuart.tau1_SP           = 1/(-Cramer.cuart.Autovalores(1));
Cramer.cuart.tau2_SP           = 1/(-Cramer.cuart.Autovalores(2));
% Cramer.cuart.w_nF              = abs(Cramer.cuart.Autovalores(3));
% Cramer.cuart.tsiF              =-real(Cramer.cuart.Autovalores(3))/Cramer.cuart.w_nF ;
% No tiene corto periodo;

%% CRAMER de cada modo y furncion de transferencia

% Elevador a Angulo de ataque

K  = MS;
KK = ML;
KKK =MSS;
K(:,1) =0 ;
KK(1,1)= C_zdeltae;
KK(2,1)= -C_mdeltae;


Cramer.ataque.sist    = simplify(s^2*KKK+s*K+KK);
Cramer.ataque.ec      = det(Cramer.ataque.sist);
Cramer.ataque.coef    = sym2poly(Cramer.ataque.ec);
Cramer.ataque.raices  = roots(Cramer.ataque.coef);
% Cramer.Elev_cabeceo.Tau1     = -1/Cramer.Elev_cabeceo.raices(1);
% Cramer.Elev_cabeceo.Tau2     = -1/Cramer.Elev_cabeceo.raices(2);
% Cramer.Elev_cabeceo.K        = Cramer.Elev_cabeceo.coef(1)*(-1/Cramer.Elev_cabeceo.Tau1)*(-1/Cramer.Elev_cabeceo.Tau2)/((-1/Cramer.cuart.tau1_SP)*(-1/Cramer.cuart.tau2_SP)*Cramer.cuart.w_nF^2*Cramer.cuart.Cuarticalong_coef(1));


TF.ataqueNum =  tf([1],Cramer.cuart.Cuarticalong_coef);
TF.ataqueDen =  tf(Cramer.ataque.coef, [1]);
TF.ataque    =  TF.ataqueDen*TF.ataqueNum;



% Elevador a Velocidad Angular de Cabeceo
K  = MS;
KK = ML;
KKK =MSS;
K(:,2) =0 ;
KK(1,2)= C_zdeltae;
KK(2,2)= -C_mdeltae; %%%%%%%%%%%%%%%Cmd deltae


Cramer.Elev_cabeceo.sist    = simplify(s^2*KKK+s*K+KK);
Cramer.Elev_cabeceo.ec      = det(Cramer.Elev_cabeceo.sist);
Cramer.Elev_cabeceo.coef    = sym2poly(Cramer.Elev_cabeceo.ec);
Cramer.Elev_cabeceo.raices  = roots(Cramer.Elev_cabeceo.coef);
% Cramer.Elev_cabeceo.Tau1     = -1/Cramer.Elev_cabeceo.raices(1);
% Cramer.Elev_cabeceo.Tau2     = -1/Cramer.Elev_cabeceo.raices(2);
% Cramer.Elev_cabeceo.K        = Cramer.Elev_cabeceo.coef(1)*(-1/Cramer.Elev_cabeceo.Tau1)*(-1/Cramer.Elev_cabeceo.Tau2)/((-1/Cramer.cuart.tau1_SP)*(-1/Cramer.cuart.tau2_SP)*Cramer.cuart.w_nF^2*Cramer.cuart.Cuarticalong_coef(1));


TF.vel_cabeceoNum =  tf([1],Cramer.cuart.Cuarticalong_coef);
TF.vel_cabeceoDen =  tf(Cramer.Elev_cabeceo.coef, [1]);
TF.vel_cabeceo    =  TF.vel_cabeceoDen*TF.vel_cabeceoNum;
TF.Asiento        =  TF.vel_cabeceo*tf([1],[1 0]);

% Funcion de transferencia del actuador Tlag = 0.05
%Primer orden

TF.actuador = tf([1],[0.05 1]);

%Segundo orden: Sacado del actuador del modelo ' aeroblk_HL20 '
wn_act=44;
z_act=   0.707106781186547;
TF.actuador = tf([1],[1/wn_act^2 2*z_act/wn_act 1]);

%% Funciones de transferencia del trabajo 1

TF_T1=load('TF_OL.mat');

% Guardamos los valores en TF

TF.Den =      tf(TF_T1.TF.long.Den,[1]);
TF.Nude =     tf(TF_T1.TF.long.Nums(1,:),[1]);
TF.Nalfade =  tf(TF_T1.TF.long.Nums(2,:),[1]);
TF.Nthetade = tf(TF_T1.TF.long.Nums(3,:),[1]);
TF.Nqde     = tf(TF_T1.TF.long.Nums(4,:),[1]);

TF.vel_cabeceoNum =  tf([1],Cramer.cuart.Cuarticalong_coef);
TF.vel_cabeceoDen =  tf(Cramer.Elev_cabeceo.coef, [1]);
TF.vel_cabeceo    =  TF.Nqde/TF.Den;
TF.Asiento        =  TF.Nthetade/TF.Den;


%% Manual Estacionario
 TimeEnd=120; %segundos
 StepSize=0.01;
 vTime = [0:StepSize:TimeEnd];
%  titles={'Respuesta temporal Velocidad vs Entrada escalón elevador','Respuesta temporal \alpha vs Entrada escalón elevador',...
%      'Respuesta temporal \theta vs Entrada escalón elevador', 'Respuesta temporal q vs Entrada escalón elevador'};
%  ylabels={'$$\hat{u}$$','$$\alpha$$ ( $$^{o}$$ )','$$\theta$$ ( $$^{o}$$ )','q (rad/s)'};
 vOut=step(TF.Asiento,vTime);
 vEstacionario=[0,0,0,0];
 overshoot_long=[0,0,0,0];
 undershoot_long=[0,0,0,0];

    figure(1)
    plot(vTime,vOut,'b')
    grid on; hold all;
    xlabel('Tiempo(s)');
%     ylabel(ylabels(i),'Interpreter','Latex');
%     title(titles(i));
    vEstacionario=vOut(length(vTime));
    plot(vTime,ones(size(vTime))*vEstacionario,'--','LineWidth',1);
    
%     overshoot_long(i)= -vEstacionario(i)+max(vOut(:,:,i));
%     undershoot_long(i)= -vEstacionario(i)+min(vOut(:,:,i));
%      
%     format_Grafico = strcat('EscalonEstacionario',num2str(i));
%     format_Grafico = [pwd filesep 'Figuras' filesep format_Grafico];
%     saveas(gcf,format_Grafico,'epsc'); 

vEstacionario=TF.Nthetade.Numerator{1}(end)/TF.Den.Numerator{1}(end);

%% Funcion transferencia Direct Link

K_DL=  1/vEstacionario;%1/vEstacionario;
TF.DirectLink = tf([K_DL], [1]);

%Funcion de transferencia global
TF.Global = TF.DirectLink*TF.actuador*TF.Asiento;

figure(2)
 vOut=step(TF.Global,vTime);
hold all
       plot(vTime,ones(size(vTime))*1,'--r','LineWidth',1);
            plot(vTime,vOut,'b'); 

 grid on; hold all;
    xlabel('Tiempo(s)');
    ylabel('\Delta\theta ','Interpreter','Tex');
    title('Respuesta temporal a escalón unitario en \Delta\delta_s_t_i_c_k');
        setFigureDefaults()
format_Grafico = 'EscalonEstacionario';
    format_Grafico = [pwd filesep 'Figuras' filesep format_Grafico];
    saveas(gcf,format_Grafico,'epsc'); 

%% Respuesta de las variables a rampa
TimeEnd     = 100;
StepSize    = 0.01;
vTime       = [0:StepSize:TimeEnd];  %Querry points
    vInpTime    = [0 1/5 3 TimeEnd];
% vInpTime    = [0 10 15 TimeEnd];
vInp        = [0 1 1 1];             %Sample pints values
vInput      = interp1(vInpTime,vInp,vTime);

% TFlongvector= [TF.ataque, TF.Asiento];
TFlongvector= [TF.DirectLink*TF.actuador*TF.Nude/TF.Den, TF.DirectLink*TF.actuador*TF.ataque,...
    TF.Global, TF.DirectLink*TF.actuador*TF.Nqde/TF.Den];
%Longitudinal estacionaria
for i=1:4
    vTime       = [0:StepSize:TimeEnd];
    vInpTime    = [0 1/5 3 TimeEnd]; %Sample points
    vInput      = interp1(vInpTime,vInp,vTime);

    vOut_long = lsim(TFlongvector(i),vInput,vTime);

    figure(i+2)
    titles={'Respuesta temporal u vs rampa elevador',...
         'Respuesta temporal \alpha vs rampa elevador',...
         'Respuesta temporal \theta vs rampa elevador',...
         'Respuesta temporal q vs rampa elevador'};
    ylabels={'\Deltau ','\Delta\alpha ','\Delta\theta','\Deltaq'};

    plot(vTime,vInput,'r');
    grid on;hold all;
    
    plot(vTime,vOut_long,'b');

    grid on;hold all;
    title(titles(i))
    xlabel('Tiempo(s)');
    ylabel(ylabels(i));
    legend({'Input','Output'})
    set(gcf,'Color',[1 1 1]);
    
%     vEstacionario_rampa(i)=vOut_long(length(vTime));
%     plot(vTime,ones(size(vTime))*vEstacionario_rampa(i),'--','LineWidth',1);
%     
%     overshoot_long(i)= -vEstacionario_rampa(i)+max(vOut_long);
%     undershoot_long(i)= -vEstacionario_rampa(i)+min(vOut_long);

    setFigureDefaults()
    format_Grafico = strcat('RampaEstacionario',num2str(i));
    format_Grafico = [pwd filesep 'Figuras' filesep format_Grafico];
    saveas(gcf,format_Grafico,'epsc'); 


end

%% Delay-time

% Calculo de la derivada
TF.Global_dot=TF.Global*tf([1, 0],[1]);

% Vectores de tiempo, input
StepSize=0.0001;
vTime       = [0:StepSize:TimeEnd];
vInpTime    = [0 1/5 3 TimeEnd]; %Sample points
vInput      = interp1(vInpTime,vInp,vTime);

%Obtencion de la respuesta
vOut_long = lsim(TF.Global,vInput,vTime);
vOut_long_dot = lsim(TF.Global_dot,vInput,vTime);

%Representacion grafica
figure(5)
setFigureDefaults()
titles={'Rise time, Delay time \theta vs rampa elevador','Interpreter','Latex'};
ylabels={'\Delta\theta '};

%Plot input
pinput=plot(vTime,vInput,'r');
grid on;hold all;

%Plot output y derivada
presponse=plot(vTime,vOut_long,'b');
presponsedot=plot(vTime,vOut_long_dot,'c');

%Maximo de la derivada
[y_dot_max, x_m_max] = max(vOut_long_dot);
m = y_dot_max;
x = vTime(x_m_max);
y =  vOut_long(x_m_max);

%Plot de la recta tangente
n =y-m*x;
x_tang=linspace(0,1,10000); 
y_tang=n+m*x_tang;
plot(x_tang, y_tang,'--k','LineWidth',1 )
 
%Buscar corte con el eje y, Time Delay
x_time=find(abs(y_tang)<0.001);
timeDelay=x_tang(x_time(round(length(x_time)/2)));
ptimedelay=plot(linspace(0,timeDelay,1000),zeros(1000),'g');

%Rise time
[y_01, x_01]= min(abs(vOut_long-0.1));
% [y_09, x_09]= min(abs(vOut_long-0.9)); % Not general
[y_09_vect, x_09_vect]=findpeaks(-abs(vOut_long-0.9)); 
x_09=x_09_vect(1);
y_09=y_09_vect(1);
riseTime=vTime(x_09)-vTime(x_01);
y_height=-0.1;
% plot(vTime,abs(vOut_long-0.9))
prisetime=plot(linspace(vTime(x_01),vTime(x_09),100),linspace(y_height,y_height,100),'m');
plot(vTime(x_09),y_height,'xm');plot(vTime(x_01),y_height,'xm');


grid on;hold all;
title(titles(1))
xlabel('Tiempo(s)');
ylabel(ylabels(1));
legend([pinput,presponse,presponsedot],'Input','Step Response','Step Response dot','Time Delay','Rise Time','location','best','Interpreter','tex')
set(gcf,'Color',[1 1 1]);

%Plot puntos relevantes
plot( timeDelay,0,'k+','LineWidth',1);
plot(vTime(x_m_max),y_dot_max,'+');
plot(linspace(vTime(x_01),vTime(x_09),100),linspace(y_height,y_height,100),'m')

plot(vTime(x_09),y_height,'xm');plot(vTime(x_01),y_height,'xm');
plot(linspace(vTime(x_01),vTime(x_01),100),linspace(y_height,0.1,100),'-k','linewidth',0.1)
plot(linspace(vTime(x_09),vTime(x_09),100),linspace(y_height,0.9,100),'-k','linewidth',0.1)
    
%     vEstacionario_rampa(i)=vOut_long(length(vTime));
%     plot(vTime,ones(size(vTime))*vEstacionario_rampa(i),'--','LineWidth',1);
%     
%     overshoot_long(i)= -vEstacionario_rampa(i)+max(vOut_long);
%     undershoot_long(i)= -vEstacionario_rampa(i)+min(vOut_long);


    xlim([0,vTime(x_09)+0.5])
ylim([-0.5,1.5])
setFigureDefaults()
format_Grafico = strcat('DelayTime');
format_Grafico = [pwd filesep 'Figuras' filesep format_Grafico];
saveas(gcf,format_Grafico,'epsc'); 

    
% %%
% F_alfa = [-0.5 -0.2 0.0 0.5 1.0 2.0 3.0 4.0];
% F_q    = [-0.5 -0.2 0.0 0.5 1.0 2.0 3.0 4.0];
% 

% C_mq_dim=C_mq*Geo.c/(2*FC.us);
% 
% kdealfa=-(F_alfa-1)*C_malfa/C_mdeltae;

% kdeq=-(F_q-1)*C_mq_dim/C_mdeltae;
% 




