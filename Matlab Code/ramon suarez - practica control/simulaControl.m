%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      %
%    Diseño de Aplicaciones basadas en Drones 2016.    %
%               Organiza Vision4UAV-UPM                % 
%      http://www.vision4uav.com/SummerSchool16        %
%                                                      %
%  Autor:                                              %
%  Ramon A. Suarez Fernandez (suarez.ramon@gmail.com)  %
%                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulaci�n controlador

%ruta=[12 0;12 5;7 -5;-1 4;-6 -5;-2 0];
%ruta=[0 0; 0 20; 0 10; 10 10; -10 10; 0 10; 2 10; 0 12; -2 10; 0 8];
ruta = [4 10; 0 14; -4 10; 0 6; 6 10; 0 16; -6 10; 0 4; 8 10; 0 18];

radio=0.05;
sim('simuladorControlCuatrirrotor_summerschool_students')

%% Representaci�n gr�fica y c�lculo de indicadores
indicadores=pintascdata(scdata,ruta,radio);

%%
%% Evaluaci�n 

clc;
disp('************ Rendimiento GPP **************')
%Preparamos el vector para evaluar rendimiento
J=[indicadores.dmedia, indicadores.dmax, indicadores.tiempoRecorrido];
taux=indicadores.tiempoMin;

%Preparamos la matriz de preferencias
%          [-----AD----](-----D----](-----T----](-----I----](-----AI---](---->
PhyMatrix=[0         0.0500      0.1000      0.1500      0.2500      0.500;
           0         0.1000      0.1500      0.2500      0.5000      1.000;
           taux         2*taux    2.5*taux      4*taux      5*taux     9*taux];

Etiquetas={'Distancia Media','Distancia M�xima','Tiempo de Recorrido';
               'metros',           'metros',        'segundos'};
       
%Evaluamos

[GPP, Porcentajes, Zonas]=Rendimiento_GPP(J,PhyMatrix,Etiquetas)
