function [ score ] = rootLocusScore_Carlos(x)
    
    global TF
    Kdealpha = x(1)*1e0;
    Kdeq     = x(2)*1e-3;
    TFCL.theta = (TF.Actuator*TF.DirectLink*TF.Nthetade)/(TF.Den+TF.Actuator*(TF.SensorAlpha*Kdealpha*TF.Nalfade+TF.Sensorq*Kdeq*TF.Nqde));
    TFCL.TFOL  =  TF.Actuator * (TF.SensorAlpha*Kdealpha*TF.Nalfade + TF.Sensorq*Kdeq*TF.Nqde) / TF.Den;

%Ordenar por la parte real de los polos
    aux     = pole(TFCL.theta);
    sorted  = esort(aux);

%Ordenar por el margen de ganancia
    warning('off','Control:analysis:MarginUnstable')
%     isStable = allmargin(TFCL.TFOL);
    [Gm,~,~,~] = margin(TFCL.TFOL);
    
    if real(sorted(1))>0
        score = 1e6*real(sorted(1)) - 20*log10(Gm);
    else
        score = - 20*log10(Gm);
    end

    
    
%     if isStable.Stable == 1
%         [Gm,~,~,~] = margin(TFCL.TFOL);
%         score = -20*log10(Gm);
%     else
%         score = 1e2;
%     end
    
    
    
%     if Pm<25
%         score = NaN;
%     elseif isinf(20*log10(Gm))
%         score = NaN;
%     else
%         score = -20*log10(Gm);
%     end
    
    warning('on','Control:analysis:MarginUnstable')

end

