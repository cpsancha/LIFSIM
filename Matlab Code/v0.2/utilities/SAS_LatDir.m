%*********************APARTADO 4**************************************
% 
% A partir del diagrama de bloques del enunciado, se han obtenido en el
% archivo Apartado_4_desarrollo.nb de Wolfram Mathematica las funciones de
% transferencia en lazo cerrado. A partir de ellos y utilizando los datos
% de las funciones de transferencia del trabajo 1, las funciones de 
% transferencia del actuador y del sensor, y Kdl obtenido en el apartado 3
% se obtienen las expresiones de la Funcion de transferencia en lazo
% cerrado. 

% clc
close all
clear all
global TF


%% Funcion de transferencia del actuador Tlag = 0.05
%Primer orden

TF.actuador = tf([1],[0.05 1]);

%Segundo orden: Sacado del actuador del modelo ' aeroblk_HL20 '
wn_act=44;
z_act=   0.707106781186547;
TF.actuador = tf([1],[1/wn_act^2 2*z_act/wn_act 1]);

%TF.actuador = tf (1,1);

%% Funcion de transferencia del sensor Tlag = 0.1
%Primer orden
Tlag = 0.1;
TF.sensor = tf([1],[Tlag 1]);
%TF.sensor = tf(1,1);

%% LIBIS Transfer functions

TF_T1=load('TF_OL.mat');
TF.Den =      tf(TF_T1.TF.lat.Den_a,[1]);


% TF_T1 = load('LibisTransferFunctions.mat');
% TF.Den =      tf(TF_T1.LibisTF.lat.Denominator,[1]);

% TF.Nude =     tf(TF_T1.LibisTF.long.Numerators(1),[1]);
% TF.Nalfade =  tf(TF_T1.LibisTF.long.Numerators(2),[1]);
% TF.Nthetade = tf(TF_T1.LibisTF.long.Numerators(3),[1]);
% TF.Nqde     = tf(TF_T1.LibisTF.long.Numerators(4),[1]);



%% Funcion transferencia Direct Link
% vEstacionario=TF.Nthetade.Numerator{1}(end)/TF.Den.Numerator{1}(end);
% vEstacionario=-4.4217;
K_DL=  1;%1/vEstacionario; No usamos el valor estacionario puesto que es inestable
TF.DirectLink = tf([K_DL], [1]);

%% Pole map initialization
figure()
    T1 = tf([1],[1 50]); 
    T2 = tf([1],[1 -50]);
    T3 = tf([1],[1 -50i]);
    T4 = tf([1],[1 50i]);
    T  = T1*T2*T3*T4;
    pzmap(T,'w');
    grid on
    hold on
    
    
    marker = ['+','o','*','s','^','p','v','<'];
    color  = ['k','m','c','r','g','b','y','m'];
    titles = [-0.5,-0.2,0,0.5,1,2,3,4];
    title('Barrido ganancias de realimentacion','interpreter','tex')

%% % Polos del sistema en Open Loop SIN AUMENTAR. (sin actuadores, etc)
 
         for k=1:length(pole(1/TF.Den))
           aux = pole(1/TF.Den); 
           pOL(k)=plot(real(aux(k)),imag(aux(k)),strcat('k','.'));hold on 
         end
