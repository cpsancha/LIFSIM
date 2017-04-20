function [ eigenStruct ] = getEigenData( eigenvalue )
%GETEIGENDATA Summary of this function goes here
%   Detailed explanation goes here

if length(eigenvalue)~=1
    eigenStruct.mode = cell(1,length(eigenvalue));
    eigenStruct.value(1:length(eigenvalue))   = NaN;
    eigenStruct.t_12(1:length(eigenvalue))    = NaN;
    eigenStruct.T2(1:length(eigenvalue))      = NaN;
    eigenStruct.Period(1:length(eigenvalue))  = NaN;
    eigenStruct.freqNat(1:length(eigenvalue)) = NaN;
    eigenStruct.Damp(1:length(eigenvalue))    = NaN;
    eigenStruct.Tau(1:length(eigenvalue))     = NaN;
end

for i=1:length(eigenvalue)
    
    if isreal(eigenvalue(i))
        eigenStruct.mode{i}  = 'Exponential';
        eigenStruct.value(i) = eigenvalue(i);
        if (eigenvalue(i) < 0)
            eigenStruct.t_12(i) = log(0.5)/eigenvalue(i);
        else
            eigenStruct.T2(i) = log(2)/eigenvalue(i);
        end
        eigenStruct.Tau(i) = -1/eigenvalue(i);
        
    else
        eigenStruct.mode{i}  = 'Oscillatory';
        eigenStruct.value(i) = eigenvalue(i);
        eigenStruct.Period(i) = 2*pi/abs(imag(eigenvalue(i)));
        if (real(eigenvalue(i)) < 0)
            eigenStruct.t_12(i) = log(0.5)/real(eigenvalue(i));
        else
            eigenStruct.T2(i) = log(2)/real(eigenvalue(i));
        end
        eigenStruct.freqNat(i) = abs(eigenvalue(i));
        eigenStruct.Damp(i) = -real(eigenvalue(i))/abs(eigenvalue(i));
    end
    
end

end

