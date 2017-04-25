function [ score ] = rootLocusScore( x)
    marker = ['+','o','*','s','^','p','v','<'];
    color  = ['k','m','c','r','g','b','y','m'];
gainValues = [0,0];
gainsEstable = [0,0];

global TF
kdealfa=x(1)/100;
kdeq=x(2)/100;

   for i = 1:length(kdealfa)
     for j = 1:length(kdeq)
         
         
     Kdealfa =  kdealfa(i);
     Kdeq    =  kdeq(j);
     
 
    TFCL.theta = (TF.actuador*TF.DirectLink*TF.Nthetade)/(TF.Den-TF.actuador*(TF.sensor*Kdealfa*TF.Nalfade+TF.sensor*Kdeq*TF.Nqde));
    TFCL.TFOL  = -TF.actuador*(TF.sensor*Kdealfa*TF.Nalfade+TF.sensor*Kdeq*TF.Nqde)/TF.Den;

%    [Gm.TFOL,Pm.TFOL,Wgm.TFOL,Wpm.TFOL] = margin(TFCL.TFOL) ;
%    Gm_dB.TFOL = 20*log10(Gm.TFOL);
   % getEigenData(pole(TFCL.theta))
     aux = pole(TFCL.theta); 
     auxReal = real(aux);
     sorted=esort(aux);
     score=real(sorted(1));

% ispositive = auxReal >0;
% 
% if sum(ispositive)==0
%     disp('Estable')
%     disp(sum(ispositive))
%     gainsEstable(end+1,:)= [Kdealfa, Kdeq];
% 
% elseif sum(ispositive)==1
%     disp(sum(ispositive))
%     gainValues(end+1,:)= [Kdealfa, Kdeq];
%     score=0;
%     
%     disp(gainValues)

       
% else
%     disp('Inestable')
%     disp(sum(ispositive))



% end

     
        
        
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
end

