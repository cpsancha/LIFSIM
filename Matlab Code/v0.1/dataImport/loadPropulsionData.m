%          DEFINITIONS:    

%          J=V/nD (advance ratio)                                                                               
%          Ct=T/(rho * n**2 * D**4) (thrust coef.)                                                              
%          Cp=P/(rho * n**3 * D**5) (power coef.)                                                               
%          Pe=Ct*J/Cp (efficiency)                                                                              
%          V          (model speed in MPH)             
%          PWR        (Hp)   
%          Torque     (In-Lbf)     
%          Thrust     (Lbf) 
% Perform unit conversions
%     Hp2Watts = 745.7;
%     PWR = Hp2Watts.*PWRHp;    
%     V = convvel(Vmph, 'mph', 'm/s'); 
%     Eta = Pe;

load('propeller17x8.mat')
    propeller.eta=reshape(Pe1,[30,13]);    
    propeller.V=reshape(V1,[30,13]);
    propeller.V=0.44704.*propeller.V; %to m/s
    propeller.RPM=rpmmat;
    
    for i=1:13
        rpmvec(i,:)=linspace(1000*i,1000*i,30);
        x=propeller.V(:,i);
        v=propeller.eta(:,i);
        xq=linspace(0,63,1000);
        Vq(i,:)=interp1(x,v,xq,'pchip',0);
%         plot(xq,Vq(i,:)); hold on
    end
    
    
    LD.propeller1.eta = Vq;
    LD.propeller1.V   = xq;
    LD.propeller1.RPM = rpmvec(:,1);
    
    
    % Propeller2 data (VERTICAL) 
    
    %Different rpm range, up to 15000 rpm
    propeller2.eta=reshape(Pe2,[30,15]);    
    propeller2.V=reshape(V3,[30,15]);
    propeller2.V=0.44704.*propeller2.V; %to m/s
    for i=1:15
        rpmvec(i,:)=linspace(1000*i,1000*i,30);
        x=propeller2.V(:,i);
        v=propeller2.eta(:,i);
        xq=linspace(0,83,1000);
        Vq(i,:)=interp1(x,v,xq,'pchip',0);
        plot(xq,Vq(i,:)); hold on
    end
    
    LD.propeller2.eta = Vq;
    LD.propeller2.V   = xq;
    LD.propeller2.RPM = rpmvec(:,1);
    
    
%     rpmGrid=reshape(rpmvec',[390,1]);
   

    
    
%         C = horzcat(A1,...,AN)
%     
%     propeller.J=reshape(J1,[30,13]);  
% 
%     propeller.Cp=reshape(Cp1,[30,13]);
%     propeller.Ct=reshape(Ct1,[30,13]);
%     propeller.Cp=reshape(Cp1,[30,13]);
%     propeller.PWR=745.7.*reshape(PWR1,[30,13]);
%     propeller.Thrust=4.448222.*reshape(eThrus1,[30,13]);  %lbf to Nenwtons
%     propeller.Torque=4.448222.*0.0254.*reshape(Torqu1,[30,13]);  %lbf*inc to N*m
%     
% %     vq = griddata(V1,rpmGrid,PWR1,82,11000)
%     
%     motor1.PWRI=PI;
%     motor1.RPMI=NI;
%     motor1.etaI=etaI;
%     
    