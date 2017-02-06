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

%% Funci�n distanciaARecta()

function d=distanciaARecta(pto,recta)
% Calcula la distancia de un punto a una recta
% pto=[xa,ya]
% recta=[m b]  -> y=m*x+b
%   Nota: si la recta es vertical poner pasando por x1 poner recta=[inf,x1]
% d=distanciaARecta(pto,recta)
%
if recta(1)==inf
    d=pto(1)-recta(2);
else
    d=abs(recta(1)*pto(1)-pto(2)+recta(2))/sqrt(recta(1)^2+1);
end
