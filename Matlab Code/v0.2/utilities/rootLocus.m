function [ TFCL ] = rootLocus( Kdealfa,Kdeq, TF )

       
    TFCL.u    = (TF.actuador*TF.DirectLink*TF.Nude)/(TF.Den-TF.actuador*(TF.sensor*Kdealfa*TF.Nalfade+TF.sensor*Kdeq*TF.Nqde));
    TFCL.alfa  = (TF.actuador*TF.DirectLink*TF.Nalfade)/(TF.Den-TF.actuador*(TF.sensor*Kdealfa*TF.Nalfade+TF.sensor*Kdeq*TF.Nqde));
    TFCL.theta = (TF.actuador*TF.DirectLink*TF.Nthetade)/(TF.Den-TF.actuador*(TF.sensor*Kdealfa*TF.Nalfade+TF.sensor*Kdeq*TF.Nqde));
    TFCL.q     = (TF.actuador*TF.DirectLink*TF.Nqde)/(TF.Den-TF.actuador*(TF.sensor*Kdealfa*TF.Nalfade+TF.sensor*Kdeq*TF.Nqde));
    TFCL.TFOL  = -TF.actuador*(TF.sensor*Kdealfa*TF.Nalfade+TF.sensor*Kdeq*TF.Nqde)/TF.Den;


   %% getEigenData(pole(TFCL.theta))

        for k=1:length(pole(TFCL.theta))
           aux = pole(TFCL.theta); 
           
           disp(
           if real(aux)>0
               disp('Pole with real positive part')
           end
           %% p1(k)=plot(real(aux(k)),imag(aux(k)),strcat(color(i),marker(j)));hold on 
        end
end

