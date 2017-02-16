%% *******************************************************************
%  *                                                                 *
%  *  Script to perform various data checks.                         *
%  *                                                                 *
%  *  Documentation to be coded...                                   *
%  *                                                                 *
%  *******************************************************************


%% DEFINE NON EXISTANCE AND FLAGS
checks.nonExist  = 9999; %Value of the added field if the array has only one component
checks.isEmptyFlag  = 1; %Value of the flag in the case of empty, must be greater than nonEmptyFlag
checks.nonEmptyFlag = 0; %Value of the flag in the case of not empty
checks.consistentLength   = 0;    %Value of the flag in the case of consistent length
checks.inconsistentLength = 1;    %Value of the flag in the case of inconsistent length, must be greater than ConsistentLength
checks.consistentDimension   = 0; %Value  of the flag in the case of consistent dimension
checks.inconsistentDimension = 1; %Value  of the flag in the case of inconsistent dimension, must be greater than ConsistentDimension
checks.errorMessage = 'Press Ctrl+C to stop the program, or any other key to continue.';




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
    eval(sprintf('LD.flags.%s.empty=checks.nonEmptyFlag',strCases));
    eval(sprintf('LD.flags.%s.length=checks.consistentLength',strCases));
    eval(sprintf('LD.flags.%s.dimension=checks.consistentDimension',strCases));
    %Check if empty
    if isempty(eval(sprintf('LD.%s',strCases)))
        eval(sprintf('LD.flags.%s.empty = checks.isEmptyFlag',strCases));
        eval(sprintf('wrn = msgbox(''LD.%s==[]'',''Warning'',''warn'');',strCases));
        uiwait(wrn);
        disp(checks.errorMessage)
        pause
    %Check if not vector
    elseif ~isvector(eval(sprintf('LD.%s',strCases)))
        eval(sprintf('LD.flags.%s.dimension=checks.inconsistentDimension',strCases));
        eval(sprintf('wrn = msgbox(''LD.%s has inconsistent dimensions. Must be scalar or vector.'',''Warning'',''warn'');',strCases));
        uiwait(wrn);
        disp(checks.errorMessage)
        pause
    %Check if length==1
    elseif isequal(length(eval(sprintf('LD.%s',strCases))),1)
        eval(sprintf('LD.flags.%s.length=checks.inconsistentLength;',strCases)) %Width==1 --> no se puede usar pre-lookup
        eval(sprintf('LD.%s(2) = checks.nonExist;',strCases));
    end
end




%% CHECK STABILITY COEFFICIENTS
for i=1:length(checks.coeffs)
    strCoeffs=checks.coeffs(i);
    %Clear flag
    eval(sprintf('clear LD.flags.%s',strCoeffs));
    %Set standard values for flags
    eval(sprintf('LD.flags.%s.empty=checks.nonEmptyFlag;',strCoeffs));
    for j=1:length(checks.analysisCases)
        strCases=checks.analysisCases(j);
        eval(sprintf('LD.flags.%s.length.%s=checks.consistentLength;',strCoeffs,strCases));
    end
    %Check if empty
    if isempty(eval(sprintf('LD.%s',strCoeffs)))
        eval(sprintf('LD.flags.%s.empty = checks.isEmptyFlag;',strCoeffs));
        eval(sprintf('wrn = msgbox(''LD.%s==[]'',''Warning'',''warn'');',strCoeffs));
        uiwait(wrn);
        disp(checks.errorMessage)
        pause
    else
        eval(sprintf('LD.flags.%s.empty=checks.nonEmptyFlag;',strCoeffs));
    end
    %Check if dimensions are consistent with the number of analysis cases
    if ~isequal(length(size(eval(sprintf('LD.%s',strCoeffs)))),length(checks.analysisCases))
        eval(sprintf('LD.flags.%s.dimension.numElements=checks.inconsistentDimension;',strCoeffs));
        eval(sprintf('wrn = msgbox(''LD.%s has more dimensions than cases of analysis'',''Warning'',''warn'');',strCoeffs));
        uiwait(wrn);
        disp(checks.errorMessage)
        pause
    else
        eval(sprintf('LD.flags.%s.dimension.numElements=checks.consistentDimension;',strCoeffs));
    end
    
    
    %% HASTA AQUI VA BIEN... MAÑANA MAS
    
    
    %Check if dimension is consistent with each analysis case
    for j=1:length(checks.analysisCases)
        strCases = checks.analysisCases(j);
        if ~isequal(size(eval(sprintf('LD.%s',strCoeffs)),j),length(eval(sprintf('LD.%s',strCases))))
            eval(sprintf('LD.flags.%s.dimension.%s=checks.inconsistentDimension;',strCoeffs,strCases));
            eval(sprintf('wrn = msgbox(''LD.%s has inconsistent dimensions with LD.%s.'',''Warning'',''warn'');',strCoeffs,strCases));
            uiwait(wrn);
            disp(checks.errorMessage)
            pause
        else
            eval(sprintf('LD.flags.%s.dimension.%s=checks.consistentDimension;',strCoeffs,strCases));
        end
    end
    %Check if length==1
    if isequal(length(eval(sprintf('LD.%s',strCoeffs))),1)
        eval(sprintf('LD.flags.%s.length=checks.inconsistentLength',strCoeffs)) %Width==1 --> no se puede usar pre-lookup
        eval(sprintf('LD.%s(2) = checks.nonExist;',strCoeffs));
    end
        
end        
        
clear i j strCases strCoeffs        
        
        
        
        
        
        
        
        
        
LD.flags.(checks.analysisCases{4});

eval(sprintf('if ~isequal(size(LD.CD0,i),length(LD.%s)),LD.flags.CD0.%s=inchecks.ConsistentLength,end',strCases,strCases))
 









