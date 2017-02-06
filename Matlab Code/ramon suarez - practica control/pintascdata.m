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


function indicadores=pintascdata(scdata,ruta,radio)
% indicadores=pintascdata(scdata,ruta,radio)
% 
% Dibuja las graficas de un ensayo de control
% y devuelve la informacion de los indicadores para la evaluacion
%
%    indicadores.ruta: matriz con las coord de la ruta a seguir
%    indicadores.longitudRuta: longitud de la ruta a seguir
%    indicadores.velMax: velocidad maxima del equipo (estimado)
%    indicadores.tiempoMin: tiempo minimo del recorrido (estimado)
%    indicadores.distancia: distancia de la trayectoria real en cada
%                   periodo de muestreo a la recta que une dos puntos 
%                   consecutivos de la ruta a seguir.
%    indicadores.dmedia: distancia media de la trayectoria real a la 
%                   recta que une los puntos de la ruta a seguir.
%    indicadores.dmax: distancia maxima de la trayectoria real a la 
%                   recta que une los puntos de la ruta a seguir.
%    indicadores.tiempoRecorrido: tiempo empleado en el recorrido real.
%
%
% scdata : datos obtenidos de la simulacion
% ruta : matriz con las coordenadas de los puntos a alcanzar
% radio : radio de la circunferencia donde se considera que se ha alcanzado
%   el punto de destino.

%% calculo auxiliar para pintar zona de aproximacion a cada punto de
% llegada.
theta=0:pi/10:2*pi;
x=radio*cos(theta);
y=radio*sin(theta);

%% A�ade puto inicial a la ruta, todas las rutas parten de [0,0]
ruta=[0,0;ruta];
indicadores.ruta=ruta;

%% Calcula longitud de la ruta
tramos=ruta(2:end,:)-ruta(1:end-1,:);
longitudRuta=0;
for i=1:size(tramos,1)
    longitudRuta=longitudRuta+norm(tramos(i,:));
end
indicadores.longitudRuta=longitudRuta;

%% Se supone velocidad maxima 3 m/sec
velMax=3;
indicadores.velMax=velMax;

%% Tiempo minimo estimado del recorrido
tiempoMin=longitudRuta/velMax;
indicadores.tiempoMin=tiempoMin;

%% calculo de las distancias a la trayectoria minima
distancia=zeros(size(scdata.time));
ptoini=ruta(1,:);
ptofin=ruta(2,:);
dx=ptofin(1)-ptoini(1);
dy=ptofin(2)-ptoini(2);
if dx==0
    m=inf;
    b=ptofin(1);
else
    m=dy/dx;
    b=(ptoini(2)*ptofin(1)-(ptofin(2)*ptoini(1)))/dx;
end


for i=1:size(scdata.time,1)
    
    p1=scdata.signals(1,1).values(i,3:4);
    pto=scdata.signals(1,1).values(i,1:2);
    if sum(p1~=ptofin)
        ptoini=ptofin;
        ptofin=p1;
        dx=ptofin(1)-ptoini(1);
        dy=ptofin(2)-ptoini(2);
        if dx==0
            m=inf;
            b=ptofin(1);
        else
            m=dy/dx;
            b=(ptoini(2)*ptofin(1)-(ptofin(2)*ptoini(1)))/dx;
        end
    end
    distancia(i)=distanciaARecta(pto,[m,b]);
end
indicadores.distancia=distancia;

%% calula distancia media y maxima
dmedia=sum(distancia)/size(distancia,1);
dmax= max(distancia);
indicadores.dmedia=dmedia;
indicadores.dmax=dmax;
indicadores.tiempoRecorrido=scdata.time(end);

%% Representaci�n gr�fica XY
figure
% dibuja los zonas de llegada a cada punto
for i=1:size(ruta,1)
    if i==size(ruta,1)
        color='g';
    else
        color='r';
    end
    hold on,fill(x+ruta(i,1),y+ruta(i,2),color);
end
% une los puntos de la ruta 
hold on,plot([0;ruta(:,1)],[0;ruta(:,2)],'r')
% dibuja la trayectoria seguida
hold on,plot(scdata.signals(1,1).values(:,1),scdata.signals(1,1).values(:,2))
xlabel('x')
ylabel('y')
title(['Tiempo del recorrido: ' num2str(scdata.time(end)) ' segundos.'])

%% Representaci�n gr�fica se�ales ensayo
figure
subplot(311)
plot(scdata.time,scdata.signals(1,1).values)
legend('x','y')
title(['Tiempo del recorrido: ' num2str(scdata.time(end)) ' segundos.'])
subplot(312)
plot(scdata.time,scdata.signals(1,2).values)
legend('Ref\_pitch','Ref\_roll')
subplot(313)
plot(scdata.time,distancia)
hold on,plot(scdata.time([1 end]),[dmedia dmedia],'g')
hold on,plot(scdata.time([1 end]),[dmax dmax],'r')
legend('dist.','dist. media','dist. max')
xlabel('tiempo (seg)')
ylabel('distancia (m)')
title(['Distancia a la trayectoria m�nima. Dist. media: ' num2str(dmedia) '. Dist. max: ' num2str(dmax)])
