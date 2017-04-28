function [ gainsEstable, gainValues, TF ] = rootLocus( kdealfa,kdeq, TF )
    marker = ['+','o','*','s','^','p','v','<'];
    color  = ['k','m','c','r','g','b','y','m'];
gainValues = [0,0];
gainsEstable = [0,0];
   for i = 1:length(kdealfa)
     for j = 1:length(kdeq)
         
         
     Kdealfa =  kdealfa(i);
     Kdeq    =  kdeq(j);
     
    TFCL.u    = (TF.actuador*TF.DirectLink*TF.Nude)/(TF.Den-TF.actuador*(TF.sensor*Kdealfa*TF.Nalfade+TF.sensor*Kdeq*TF.Nqde));
    TFCL.alfa  = (TF.actuador*TF.DirectLink*TF.Nalfade)/(TF.Den-TF.actuador*(TF.sensor*Kdealfa*TF.Nalfade+TF.sensor*Kdeq*TF.Nqde));
    TFCL.theta = (TF.actuador*TF.DirectLink*TF.Nthetade)/(TF.Den-TF.actuador*(TF.sensor*Kdealfa*TF.Nalfade+TF.sensor*Kdeq*TF.Nqde));
    TFCL.q     = (TF.actuador*TF.DirectLink*TF.Nqde)/(TF.Den-TF.actuador*(TF.sensor*Kdealfa*TF.Nalfade+TF.sensor*Kdeq*TF.Nqde));
    TFCL.TFOL  = -TF.actuador*(TF.sensor*Kdealfa*TF.Nalfade+TF.sensor*Kdeq*TF.Nqde)/TF.Den;


   % getEigenData(pole(TFCL.theta))
aux = pole(TFCL.theta); 
auxReal = real(aux);

ispositive = auxReal >0;

if sum(ispositive)==0
    disp('Estable')
    disp(sum(ispositive))
    gainsEstable(end+1,:)= [Kdealfa, Kdeq];
elseif sum(ispositive)==1
    disp(sum(ispositive))
    gainValues(end+1,:)= [Kdealfa, Kdeq];
    
%     disp(gainValues)

       
% else
%     disp('Inestable')
%     disp(sum(ispositive))



end

        for k=1:length(pole(TFCL.theta))
           aux = pole(TFCL.theta); 
           p1(i,j,k)=plot(real(aux(k)),imag(aux(k)),strcat(color(i),marker(j)));hold on 
        end
        
        
%         for k=1:length(pole(TFCL.theta))
%            
% 
% 
%            
%            
% 
%            if real(aux(k))>0
%                disp('Pole with real positive part')
%                disp( real(aux(k)))
%            end
%            %% p1(k)=plot(real(aux(k)),imag(aux(k)),strcat(color(i),marker(j)));hold on 
%         end
        
      end
   end
   plot(gainsEstable(:,1),gainsEstable(:,2),'+')
end

