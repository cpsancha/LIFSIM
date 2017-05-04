function [ coeff ] = fillStabilityCoeffs( LD, dataCases, Data ) %#ok<STOUT>
%FILLSTABILITYCOEFFS Summary of this function goes here
%   Detailed explanation goes here
    
    %Get indexes of the analysis data cases
    indexes = [];
    sizes   = [];
    for i=1:length(dataCases)
        index = find(contains(LD.Stability.analysisCases,dataCases{i}));
        indexes(end+1) = index; %#ok<AGROW>
        sizes(end+1) = length(LD.(LD.Stability.analysisCases{index})); %#ok<AGROW>
    end
    
    %Check if vector for reshaping
    if isvector(Data) && isequal(length(dataCases),1)
        Data = reshape(Data,[length(Data), 1]);
        sizes(end+1) = 1;
    end
    
    %Check consistent dimensions
    if ~isequal(size(squeeze(Data)),sizes)
        wrn = msgbox({'The size of the provided derivative data is inconsistent with the dimensions of the analysis cases.',...
                      'The program should be closed and fixed the derivative data or unexpected results could be obtained.'},'Warning','warn');
        uiwait(wrn);
        disp('Press Ctrl+C to stop the program, or any other key to continue.')
        pause
    end
    
    
    %Prepare pattern
    pattern = '(';
    for j=1:length(LD.Stability.analysisCases)
        if sum(contains(string(indexes),string(j)))
            pattern = strcat(pattern,':');
        else
            pattern = strcat(pattern,'1');
        end
        if ~isequal(j,length(LD.Stability.analysisCases))
            pattern = strcat(pattern,',');
        end
    end
    pattern = strcat(pattern,')');
    
    eval(sprintf('coeff%s=Data;',pattern));
    
end

