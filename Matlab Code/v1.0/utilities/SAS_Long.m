%% CLEAR VARIABLES
close all
clearvars -except LD optimresults
clc
global TF

%% ACTUATOR TRANSFER FUNCTION
%Primer orden --> Tlag = 0.05
    TF.Actuator = tf(1,[LD.LFO.timeLag 1]);

%Segundo orden: Sacado del actuador del modelo ' aeroblk_HL20 '
    wn_act = LD.Actuators.LSO.naturalFreq;
    z_act  = LD.Actuators.LSO.dampingRatio;
    TF.Actuator = tf(1,[1/wn_act^2 2*z_act/wn_act 1]);
    clear wn_act z_act

    

%% SENSOR TRANSFER FUNCTION
%Primer orden --> Tlag = 0.1
    Tlag = 0.01;
    TF.SensorAlpha = tf(1,[Tlag 1]);
    TF.Sensorq     = tf(1,[Tlag 1]);
    clear Tlag
    


%% LIBIS TRANSFER FUNCTIONS

load('LibisTransferFunctions.mat');

TF.Den =      tf(LibisTF.long.Denominator,1);
TF.Nude =     tf(LibisTF.long.Numerators(1),1);
TF.Nalfade =  tf(LibisTF.long.Numerators(2),1);
TF.Nthetade = tf(LibisTF.long.Numerators(3),1);
TF.Nqde     = tf(LibisTF.long.Numerators(4),1);



%% DIRECT LINK TRANSFER FUNCTION
    vEstacionario = LibisTF.long.StaticGains.Elevator_Theta; %#ok<NASGU>
    K_DL = 1;%1/vEstacionario; No usamos el valor estacionario puesto que es inestable
    TF.DirectLink = tf(K_DL,1);
    clear vEstacionario K_DL
    
    

%% CLOSED-LOOP TRANSFER FUNCTIONS
% TF-CL definiendo valores de Ganancias de realimentacion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Barrido de ganancias estables:
    KdealphaVect =  0.240807761235244*1e+0;
    KdeqVect     =  0.614669052712092*1e-4;
%     KdealphaVect = optimresults.x(1)*1e+0;
%     KdeqVect     = optimresults.x(2)*1e+0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % GA Results:
%     Kdealpha =  0.567/100
%              =  0.240807761235244*1e+0;
%     Kdeq     =  1.602/100;
%              =  0.614669052712092*1e-4;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:length(KdealphaVect)
    Kdealpha = KdealphaVect(i);
    
    for j = 1:length(KdeqVect)
        Kdeq = KdeqVect(j);   

        TFCL.u(i,j)     = (TF.Actuator*TF.DirectLink*TF.Nude)/(TF.Den+TF.Actuator*(TF.SensorAlpha*Kdealpha*TF.Nalfade+TF.Sensorq*Kdeq*TF.Nqde));
        TFCL.alpha(i,j) = (TF.Actuator*TF.DirectLink*TF.Nalfade)/(TF.Den+TF.Actuator*(TF.SensorAlpha*Kdealpha*TF.Nalfade+TF.Sensorq*Kdeq*TF.Nqde));
        TFCL.theta(i,j) = (TF.Actuator*TF.DirectLink*TF.Nthetade)/(TF.Den+TF.Actuator*(TF.SensorAlpha*Kdealpha*TF.Nalfade+TF.Sensorq*Kdeq*TF.Nqde));
        TFCL.q(i,j)     = (TF.Actuator*TF.DirectLink*TF.Nqde)/(TF.Den+TF.Actuator*(TF.SensorAlpha*Kdealpha*TF.Nalfade+TF.Sensorq*Kdeq*TF.Nqde));
        TFCL.TFOL(i,j)  =  TF.Actuator * (TF.SensorAlpha*Kdealpha*TF.Nalfade + TF.Sensorq*Kdeq*TF.Nqde) / TF.Den;
    end
end

%% BARRIDO EN GANANCIAS
% Dibujamos los polos de las funciones de transferencia
    figure()

%Los puntos T1, T2, T3, T4 corresponden a los límites de los ejes a representar
    warning('off','Control:ltiobject:TFComplex')
    T1 = tf(1,[1 50]); 
    T2 = tf(1,[1 -50]);
    T3 = tf(1,[1 -50i]);
    T4 = tf(1,[1 50i]);
    T  = T1*T2*T3*T4;
    clear T1 T2 T3 T4
    warning('on','Control:ltiobject:TFComplex')
    pzmap(T,'w');
    grid on
    hold on
    
%Definimos opciones graficas
    marker = ['+','o','*','s','^','p','v','<'];
    color  = ['k','m','c','r','g','b','y','m'];
    titles = [-0.5,-0.2,0,0.5,1,2,3,4];
    title('Barrido ganancias de realimentacion','interpreter','tex')
    
%Polos del sistema Closed-Loop aumentado (Incluye actuadores etc) para diferentes combinaciones de ganacias
    for i = 1:length(KdealphaVect)
        for j = 1:length(KdeqVect)
            %se pintan todos los polos para esa funcion de transferencia misma color y marcador
            for k=1:length(pole(TFCL.theta(i,j)))
                aux = pole(TFCL.theta(i,j));
                p1(i,j,k) = plot(real(aux(k)),imag(aux(k)),strcat(color(i),marker(j)));hold on  %#ok<SAGROW>
            end
        end
    end
 

%Polos del sistema en Open Loop SIN AUMENTAR. (sin actuadores, etc)
    for k=1:length(pole(1/TF.Den))
        aux = pole(1/TF.Den);
        pOL(k)=plot(real(aux(k)),imag(aux(k)),strcat('k','.'));hold on  %#ok<SAGROW>
    end


%Construcción de la leyenda:
    for j = 1:length(KdealphaVect)
        legendAlfa{j}= strcat('k_\delta_e_\alpha =  ',num2str(KdealphaVect(j))); %#ok<SAGROW>
    end

    for j = 1:length(KdeqVect)
        legendQ{j}= strcat('k_\delta_e_q =  ',num2str(KdeqVect(j))); %#ok<SAGROW>
    end

    legendVec = cat(1,'Open Loop sin aumentar',legendAlfa', legendQ');
    plotsVec       = cat(1,pOL(1),p1(:,1,1),p1(1,:,2)');

    legend(plotsVec, legendVec)
    hold off
         
%  legend([p1(1,1,1) p1(2,2,1) p1(3,3,1) p1(4,4,1) p1(5,5,1) p1(6,6,1) p1(7,7,1) p1(8,8,1)],...
%      'k_\delta_e_\alpha = 0.4299 k_\delta_e_q = 0.3027','k_\delta_e_\alpha  = 0.3439 k_\delta_e_q = 0.2422',...
%      'k_\delta_e_\alpha  = 0.2866 k_\delta_e_q = 0.2018','k_\delta_e_\alpha  = 0.1433 k_\delta_e_q = 0.1009',...
%      'k_\delta_e_\alpha  = 0.0 k_\delta_e_q = 0.0','k_\delta_e_\alpha  = -0.2866 k_\delta_e_q = -0.2018',...
%      'k_\delta_e_\alpha = -0.5731 k_\delta_e_q = -0.4036','k_\delta_e_\alpha  = -0.8597 k_\delta_e_q = -0.6055',...
%      'interpreter','tex','Location','northeast')

%  legend([p1(1,1,1) p1(2,2,1) p1(3,3,1) p1(4,4,1) p1(5,5,1) p1(6,6,1) p1(7,7,1) p1(8,8,1)],...
%      'k_\delta_e_\alpha = 0.0287 k_\delta_e_q = 0.3027','k_\delta_e_\alpha  = 0.0143 k_\delta_e_q = 0.2422',...
%      'k_\delta_e_\alpha  = 0.0086 k_\delta_e_q = 0.2018','k_\delta_e_\alpha  = 0.0057 k_\delta_e_q = 0.1009',...
%      'k_\delta_e_\alpha  = 0.0 k_\delta_e_q = 0.0','k_\delta_e_\alpha  = -0.0287 k_\delta_e_q = -0.2018',...
%      'k_\delta_e_\alpha = -0.0573 k_\delta_e_q = -0.4036','k_\delta_e_\alpha  = -0.1146 k_\delta_e_q = -0.6055',...
%      'interpreter','tex','Location','northeast')

