%% *******************************************************************
%  *                                                                 *
%  *  Script to perform various data checks.                         *
%  *                                                                 *
%  *  Documentation to be coded...                                   *
%  *                                                                 *
%  *******************************************************************


%% DEFINE NON EXISTANCE AND FLAGS
checks.nonExist  = 9999; %Value of the added field if the array has only one component
checks.existenceFlag   = 0; %Value of the flag in the case of existance (everything ok)
checks.nonExistenceFlag = 1; %Value of the flag in the case of not existence, must be greater than existenceFlag
checks.isEmptyFlag  = 1; %Value of the flag in the case of empty, must be greater than nonEmptyFlag
checks.nonEmptyFlag = 0; %Value of the flag in the case of not empty
checks.consistentLength   = 0;    %Value of the flag in the case of consistent length
checks.inconsistentLength = 1;    %Value of the flag in the case of inconsistent length, must be greater than ConsistentLength
checks.consistentDimension   = 0; %Value  of the flag in the case of consistent dimension
checks.inconsistentDimension = 1; %Value  of the flag in the case of inconsistent dimension, must be greater than ConsistentDimension
checks.errorMessage = 'Press Ctrl+C to stop the program, or any other key to continue.';
checks.closeMessage = 'Some terms are not defined, the program must close.';




%% CASES OF ANALYSIS AND COEFFICIENTS DEFINITION
checks.analysisCases = string({'alpha','alt','xcg','deltae','deltar','deltafr','deltafl'});
checks.coeffs = string({'CD0','CDalpha','CDalpha_dot','CDq','CDdeltae','CDdeltafr','CDdeltafl',...
                        'CL0','CLalpha','CLalpha_dot','CLq','CLdeltae','CLdeltafr','CLdeltafl',...
                        'Cm0','Cmalpha','Cmalpha_dot','Cmq','Cmdeltae','Cmdeltafr','Cmdeltafl',...
                        'CYbeta','CYp','CYr','CYdeltar','CYdeltafr','CYdeltafl',...
                        'Clbeta','Clp','Clr','Cldeltar','Cldeltafr','Cldeltafl',...
                        'Cnbeta','Cnp','Cnr','Cndeltar','Cndeltafr','Cndeltafl'});


%% CHECK CASES OF ANALYSIS
for i=1:length(checks.analysisCases)
    strCases=checks.analysisCases(i);
    
    %Clear flag
    eval(sprintf('clear LD.flags.%s',strCases));
    
    %Set standard values for flags
    LD.flags.(checks.analysisCases{i}).existence = checks.existenceFlag;
    LD.flags.(checks.analysisCases{i}).empty     = checks.nonEmptyFlag;
    LD.flags.(checks.analysisCases{i}).length    = checks.consistentLength;
    LD.flags.(checks.analysisCases{i}).dimension = checks.consistentDimension;
    
    %Check if exists
    if ~isfield(LD, char(checks.analysisCases(i)))
        LD.flags.(checks.analysisCases{i}).existence = checks.nonExistenceFlag;
        eval(sprintf('wrn = msgbox(''LD.%s does not exist. The program must close.'',''Error'',''error'');',strCases));
        uiwait(wrn);
        disp(checks.closeMessage)
        error(checks.closeMessage)
       
    %Check if empty
    elseif isempty(eval(sprintf('LD.%s',strCases)))
        LD.flags.(checks.analysisCases{i}).empty = checks.isEmptyFlag;
        eval(sprintf('wrn = msgbox(''LD.%s==[]'',''Warning'',''warn'');',strCases));
        uiwait(wrn);
        disp(checks.errorMessage)
        pause
        
    %Check if not vector
    elseif ~isvector(eval(sprintf('LD.%s',strCases)))
       LD.flags.(checks.analysisCases{i}).dimension=checks.inconsistentDimension;
        eval(sprintf('wrn = msgbox(''LD.%s has inconsistent dimensions. Must be scalar or vector.'',''Warning'',''warn'');',strCases));
        uiwait(wrn);
        disp(checks.errorMessage)
        pause    
    
    %Check if length==1
    elseif isequal(length(eval(sprintf('LD.%s',strCases))),1)
        LD.flags.(checks.analysisCases{i}).length=checks.inconsistentLength; %Width==1 --> no se puede usar pre-lookup
        LD.(checks.analysisCases{i})(end+1) = checks.nonExist;
    end
end

% %Check how many trailing inconsistent length flags exist
% checks.trailingUnidimensional = 0;
% for i=1:length(check.analysisCases)
%     checks.trailingUnidimensional = 
% end



%% CHECK STABILITY COEFFICIENTS
for i=1:length(checks.coeffs)
    strCoeffs=checks.coeffs(i);
    
    %Clear flag
    eval(sprintf('clear LD.flags.%s',strCoeffs));

    %Check if exists
    if ~isfield(LD.Stability, char(checks.coeffs(i)))
        LD.flags.(checks.coeffs{i}).existence = checks.nonExistenceFlag;
        eval(sprintf('wrn = msgbox(''LD.%s does not exist. The program must close.'',''Error'',''error'');',strCoeffs));
        uiwait(wrn);
        disp(checks.closeMessage)
        error(checks.closeMessage)
    else
        LD.flags.(checks.coeffs{i}).existence = checks.existenceFlag;
    end
    
    %Check if empty
    if isempty(LD.Stability.(checks.coeffs{i}))
        LD.flags.(checks.coeffs{i}).empty = checks.isEmptyFlag;
        eval(sprintf('wrn = msgbox(''LD.%s==[]'',''Warning'',''warn'');',strCoeffs));
        uiwait(wrn);
        disp(checks.errorMessage)
        pause
    else
        LD.flags.(checks.coeffs{i}).empty=checks.nonEmptyFlag;
    end
    
    %Check if dimensions are consistent with the number of analysis cases
%     if isequal(length(size(eval(sprintf('LD.%s',strCoeffs)))),length(checks.analysisCases))...
%         LD.flags.(checks.coeffs{i}).dimension.numElements=checks.consistentDimension;
%             
%     else
%         LD.flags.(checks.coeffs{i}).dimension.numElements=checks.inconsistentDimension;
%         eval(sprintf('wrn = msgbox(''LD.%s has more dimensions than cases of analysis'',''Warning'',''warn'');',strCoeffs));
%         uiwait(wrn);
%         disp(checks.errorMessage)
%         pause
%     end
    
    %Check if each dimension is consistent with each analysis case
    for j=1:length(checks.analysisCases)
        strCases = checks.analysisCases(j);
        if isequal(LD.flags.(checks.analysisCases{j}).length,checks.consistentLength)
            if ~isequal(size(LD.Stability.(checks.coeffs{i}),j),length(LD.(checks.analysisCases{j})))
                LD.flags.(checks.coeffs{i}).dimension.(checks.analysisCases{j})=checks.inconsistentDimension;
                eval(sprintf('wrn = msgbox(''LD.%s has inconsistent dimensions with LD.%s.'',''Warning'',''warn'');',strCoeffs,strCases));
                uiwait(wrn);
                disp(checks.errorMessage)
                pause
            else
                LD.flags.(checks.coeffs{i}).dimension.(checks.analysisCases{j})=checks.consistentDimension;
            end
        elseif isequal(LD.flags.(checks.analysisCases{j}).length,checks.inconsistentLength)
            if ~isequal(size(LD.Stability.(checks.coeffs{i}),j),1)
                LD.flags.(checks.coeffs{i}).dimension.(checks.analysisCases{j})=checks.inconsistentDimension;
                eval(sprintf('wrn = msgbox(''LD.%s has inconsistent dimensions with LD.%s.'',''Warning'',''warn'');',strCoeffs,strCases));
                uiwait(wrn);
                disp(checks.errorMessage)
                pause
            else
                LD.flags.(checks.coeffs{i}).dimension.(checks.analysisCases{j})=checks.consistentDimension;
            end
        else
            error('Some strange shit is happening, please debug me...')
        end
    end
    
    %Check if length==1
    for j=1:length(checks.analysisCases)
        if isequal(size(LD.Stability.(checks.coeffs{i}),1),1)
            LD.flags.(checks.coeffs{i}).length.(checks.analysisCases{j})=checks.inconsistentLength;
            LD.Stability.(checks.coeffs{i})(end+1) = checks.nonExist;
        else
            LD.flags.(checks.coeffs{i}).length.(checks.analysisCases{j})=checks.consistentLength;
        end
    end    
end        
        
%% RE-COUNT OF ERRORS FROM STORED FLAGS
%Initialize summation variables
checks.existenceErrorSum = 0;
checks.emptyErrorSum     = 0;
checks.dimensionErrorSum = 0;

%Check analysis cases
for i=1:length(checks.analysisCases)
    if isequal(LD.flags.(checks.analysisCases{i}).existence,checks.nonExist) 
        checks.existenceErrorSum = checks.existenceErrorSum + 1;
    end
    if isequal(LD.flags.(checks.analysisCases{i}).empty,checks.isEmptyFlag) 
        checks.emptyErrorSum = checks.emptyErrorSum + 1;
    end
    if isequal(LD.flags.(checks.analysisCases{i}).dimension,checks.inconsistentDimension) 
        checks.dimensionErrorSum = checks.dimensionErrorSum + 1;
    end
end

%Check stability derivatives
for i=1:length(checks.coeffs)
    if isequal(LD.flags.(checks.coeffs{j}).existence,checks.nonExist) 
        checks.existenceErrorSum = checks.existenceErrorSum + 1;
    end
    if isequal(LD.flags.(checks.coeffs{j}).empty,checks.isEmptyFlag) 
        checks.emptyErrorSum = checks.emptyErrorSum + 1;
    end
    for j=1:length(checks.analysisCases)
        if isequal(LD.flags.(checks.coeffs{j}).dimension.(checks.analysisCases{j}),checks.inconsistentDimension) 
            checks.dimensionErrorSum = checks.dimensionErrorSum + 1;
        end
    end
end

%Display warnings
if ~isequal(checks.existenceErrorSum+checks.emptyErrorSum+checks.dimensionErrorSum,0)
    eval(sprintf('wrn = msgbox({''Total Existence Errors: %i.'',''Eotal Empty Errors: %i.'',''Total Inconsistent Dimensions Errors: %i.''},''Warning'',''warn'');',...
        checks.existenceErrorSum,checks.emptyErrorSum,checks.dimensionErrorSum));
    uiwait(wrn);
end

%% CLEAR USED VARIABLES
clear i j strCases strCoeffs        
        
        






