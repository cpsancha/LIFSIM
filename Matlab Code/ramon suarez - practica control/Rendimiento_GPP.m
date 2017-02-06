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

%% Rendimiento_GPP

function [Rendimiento R Zona]=Rendimiento_GPP(J,Phy,Etiquetas)

if nargin==2
    Etiquetas={' ',' ',' '; ' ', ' ', ' '};
end

Population=size(J,1);            % N�mero de vectores objetivo.
Nobj=size(J,2);                  % N�mero de objetivos.
F=zeros(Population,Nobj);        % Inicializaci�n de F.
R=zeros(Population,Nobj);        % Inicializaci�n de R.
Zona=cell(Population,Nobj);
Rendimiento=zeros(Population,1); % Inicializaci�n de Rendimiento.
N = 1;                           % Curvatura en GPP.
                                 % N=1 --> Rectas.
                                 % N=2 --> Par�bolas....

% El offset ha sdo definido en la siguiente forma para reforzar la
% interpretabilidad del indicador de rendimiento GPP. Un Rendimiento>100
% significa que alguno de los indicadores de desempe�o se encuentra en la
% zona "Indeseable". Un Rendimiento=100, significa que todos los
% indicadores se encuentran en la frontera "Tolerable|Indeseable" para
% todos sus indicadores. Un Rendimiento=0, significa un rendimiento
% perfecto, i.e. seguimiento perfecto de la trayectoria de control en el
% tiempo de recorrido m�nimo, con un modelo que representa perfectamente al
% cuadrirrotor.

%       [-----HD----](-----D----](-----T----](-----U----](-----HU---](---->
offset=[0,       1/Nobj^2,   1/Nobj^0,     Nobj^2,     Nobj^4,    Nobj^6];
offset=100*offset/(Nobj^2*(Nobj+1));


% As� mismo, para complementar la regla heur�stica "One-Vs-Others", se ha
% modificado el offset que se suma en cada uno de los rangos de
% preferencia.

for population=1:Population
    for nobj=1:Nobj
        if J(population,nobj)>Phy(nobj,1) && J(population,nobj)<=Phy(nobj,2)     % Tenemos un valor altamente deseable
            F(population,nobj)=offset(1)*(Nobj+1)+ (offset(2)-offset(1))*((J(population,nobj)-Phy(nobj,1))/(Phy(nobj,2)-Phy(nobj,1)))^N;
            R(population,nobj)=100*(J(population,nobj)-Phy(nobj,1))/(Phy(nobj,2)-Phy(nobj,1));
            Zona{population,nobj}='AD - Altamente Deseable';
        elseif J(population,nobj)>Phy(nobj,2) && J(population,nobj)<=Phy(nobj,3) % Tenemos un valor deseable
            F(population,nobj)=offset(2)*(Nobj+1)+ (offset(3)-offset(2))*((J(population,nobj)-Phy(nobj,2))/(Phy(nobj,3)-Phy(nobj,2)))^N;
            R(population,nobj)=100*(J(population,nobj)-Phy(nobj,2))/(Phy(nobj,3)-Phy(nobj,2));
            Zona{population,nobj}='D - Deseable';
        elseif J(population,nobj)>Phy(nobj,3) && J(population,nobj)<=Phy(nobj,4) % Tenemos un valor tolerable
            F(population,nobj)=offset(3)*(Nobj+1)+ (offset(4)-offset(3))*((J(population,nobj)-Phy(nobj,3))/(Phy(nobj,4)-Phy(nobj,3)))^N;
            R(population,nobj)=100*(J(population,nobj)-Phy(nobj,3))/(Phy(nobj,4)-Phy(nobj,3));
            Zona{population,nobj}='T - Tolerable';
        elseif J(population,nobj)>Phy(nobj,4) && J(population,nobj)<=Phy(nobj,5) % Tenemos un valor indeseable
            F(population,nobj)=offset(4)*(Nobj+1)+ (offset(5)-offset(4))*((J(population,nobj)-Phy(nobj,4))/(Phy(nobj,5)-Phy(nobj,4)))^N;
            R(population,nobj)=100*(J(population,nobj)-Phy(nobj,4))/(Phy(nobj,5)-Phy(nobj,4));
            Zona{population,nobj}='I - Indeseable';
        elseif J(population,nobj)>Phy(nobj,5) && J(population,nobj)<=Phy(nobj,6) % Tenemos un valor altamente indeseable.
            F(population,nobj)=offset(5)*(Nobj+1)+ (offset(6)-offset(5))*((J(population,nobj)-Phy(nobj,5))/(Phy(nobj,6)-Phy(nobj,5)))^N;
            R(population,nobj)=100*(J(population,nobj)-Phy(nobj,5))/(Phy(nobj,6)-Phy(nobj,5));
            Zona{population,nobj}='AI - Altamente Indeseable';
        elseif J(population,nobj)>Phy(nobj,6)                                    % Tenemos un valor fuera de toda consideraci�n.
            F(population,nobj)=offset(6)*(Nobj+1);
            R(population,nobj)=100*(J(population,nobj)-Phy(nobj,5))/(Phy(nobj,6)-Phy(nobj,5));
            Zona{population,nobj}='FUERA DE TODA CONSIDERACION!!!';
        end        
    end
    Rendimiento(population,1)=sum(F(population,:));       % CONTROL 
end

%%
%% Desplegamos en pantalla:

for population=1:Population
    disp(['Vector objetivo ' num2str(population) ' : valor GPP : ' num2str(Rendimiento(population,1))]);
    Jind=J(population,:);
    Rind=R(population,:);
    for nobj=1:Nobj
       if Jind(1,nobj)<=Phy(nobj,2)
           disp(['---J' num2str(nobj) ': ' Etiquetas{1,nobj} ' = ' num2str(Jind(1,nobj)) ' ' Etiquetas{2,nobj} ' ---> AD - Altamente Deseable (' num2str(Rind(1,nobj)) '%)']);
       elseif Jind(1,nobj)>Phy(nobj,2) && Jind(1,nobj)<=Phy(nobj,3)
           disp(['---J' num2str(nobj) ': ' Etiquetas{1,nobj} ' = ' num2str(Jind(1,nobj)) ' ' Etiquetas{2,nobj} ' ---> D - Deseable (' num2str(Rind(1,nobj)) '%)']);
       elseif Jind(1,nobj)>Phy(nobj,3) && Jind(1,nobj)<=Phy(nobj,4)
           disp(['---J' num2str(nobj) ': ' Etiquetas{1,nobj} ' = ' num2str(Jind(1,nobj)) ' ' Etiquetas{2,nobj} ' ---> T - Tolerable (' num2str(Rind(1,nobj)) '%)']);
       elseif Jind(1,nobj)>Phy(nobj,4) && Jind(1,nobj)<=Phy(nobj,5)
           disp(['---J' num2str(nobj) ': ' Etiquetas{1,nobj} ' = ' num2str(Jind(1,nobj)) ' ' Etiquetas{2,nobj} ' ---> I - Indeseable (' num2str(Rind(1,nobj)) '%)']);
       elseif Jind(1,nobj)>Phy(nobj,5)
           disp(['---J' num2str(nobj) ': ' Etiquetas{1,nobj} ' = ' num2str(Jind(1,nobj)) ' ' Etiquetas{2,nobj} ' ---> AI - Altamente Indeseable (' num2str(Rind(1,nobj)) '%)']);
       end           
    end
end