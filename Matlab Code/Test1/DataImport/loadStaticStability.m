function [ S ] = loadStaticStability(Value,alpha,alphaVect,mach,machVect,...
                 alt,altVect,build,buildVect,grndht,grndhtVect,delta,deltaVect)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Define the length of the vectors and its variation
    %ALPHA
    if(isempty(alphaVect))
        alphaLength = 1;
        alphaVar    = 1;
    else
        alphaLength = length(alphaVect);
        alphaVar    = alphaVect./alpha;
    end
    %MACH
    if(isempty(machVect))
        machLength = 1;
        machVar    = 1;
    else
        machLength = length(machVect);
        machVar    = machVect./mach;
    end
    %ALT
    if(isempty(altVect))
        altLength = 1;
        altVar    = 1;
    else
        altLength = length(altVect);
        altVar    = altVect./alt;
    end
    %BUID
    if(isempty(buildVect))
        buildLength = 1;
        buildVar    = 1;
    else
        buildLength = length(buildVect);
        buildVar    = buildVect./build;
    end
    %GRNDHT
    if(isempty(grndhtVect))
        grndhtLength = 1;
        grndhtVar    = 1;
    else
        grndhtLength = length(grndhtVect);
        grndhtVar    = grndhtVect./grndht;
    end
    %DELTA
    if(isempty(deltaVect))
        deltaLength = 1;
        deltaVar    = 1;
    else
        deltaLength = length(deltaVect);
        deltaVar    = deltaVect./delta;
    end

 
%Initialize output Static Stability matrix
    S = zeros(alphaLength,...
              machLength,...
              altLength,...
              buildLength,...
              grndhtLength,...
              deltaLength);

%Build the Static Stability matrix           
    for i=1:alphaLength
        for j=1:machLength
            for k=1:altLength
                for l=1:buildLength
                    for m=1:grndhtLength
                        for n=1:deltaLength
                            S(i,j,k,l,m,n) = Value.*...
                                             alphaVar(i).*...
                                             machVar(j).*...
                                             altVar(k).*...
                                             buildVar(l).*...
                                             grndhtVar(m).*...
                                             deltaVar(n);
                        end
                    end
                end
            end
        end
    end


%Reduce the size of the output matrix removing unuseful information    
    S = squeeze(S);


% if strcmp(str,'cd')
% elseif strcmp(str,'cl')
% elseif strcmp(str,'cm')
% elseif strcmp(str,'cn')
% elseif strcmp(str,'ca')
% elseif strcmp(str,'xcp')
% elseif strcmp(str,'cla')
% elseif strcmp(str,'cma')
% elseif strcmp(str,'cyb')
% elseif strcmp(str,'cnb')
% elseif strcmp(str,'clb')
% elseif strcmp(str,'qqinf')
% elseif strcmp(str,'eps')
% elseif strcmp(str,'depsdalp')
% end
    
    
end

