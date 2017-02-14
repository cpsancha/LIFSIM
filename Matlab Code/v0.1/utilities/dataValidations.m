%% *******************************************************************
%  *                                                                 *
%  *  Script to perform various data checks.                         *
%  *                                                                 *
%  *  Documentation to be coded...                                   *
%  *                                                                 *
%  *******************************************************************


%% DEFINE NON EXISTANCE DATA
nonExist = 9999;


%% CHECK ALPHA
clear LD.flags.alpha
if isempty(LD.alpha)
    wrn = msgbox('Vector de angulos de ataque vacios: LD.alpha==[] ', 'Aviso','warn');uiwait(wrn);
    disp('Presiona Ctrl+C para detener el programa, o cualquier otra tecla para continuar con su ejecución.')
    pause
elseif length(LD.alpha)==1
    LD.alpha(2)    = nonExist;
    LD.flags.alpha.dim = 0; %Width==1, no se puede usar pre-lookup
else
    LD.flags.alpha.dim = 1;
end


%% CHECK ALTITUDE
clear LD.flags.alt
if isempty(LD.alt)
    wrn = msgbox('Vector de angulos de ataque vacios: LD.alt==[] ', 'Aviso','warn');uiwait(wrn);
    disp('Presiona Ctrl+C para detener el programa, o cualquier otra tecla para continuar con su ejecución.')
    pause
elseif length(LD.alt)==1
    LD.alt(2) = nonExist;
    LD.flags.alt.dim = 0; %Width==1, no se puede usar pre-lookup
else
    LD.flags.alt.dim = 1;
end


%% CHECK CENTER OF GRAVITY
clear LD.flags.xcg
if isempty(LD.xcg)
    wrn = msgbox('Vector de angulos de ataque vacios: LD.xcg==[] ', 'Aviso','warn');uiwait(wrn);
    disp('Presiona Ctrl+C para detener el programa, o cualquier otra tecla para continuar con su ejecución.')
    pause
elseif length(LD.xcg)==1
    LD.xcg(2) = nonExist;
    LD.flags.xcg.dim = 0; %Width==1, no se puede usar pre-lookup
else
    LD.flags.xcg.dim = 1;
end


%% CHECK DELTA OF ELEVATOR
clear LD.flags.deltae
if isempty(LD.deltae)
    wrn = msgbox('Vector de angulos de ataque vacios: LD.deltae==[] ', 'Aviso','warn');uiwait(wrn);
    disp('Presiona Ctrl+C para detener el programa, o cualquier otra tecla para continuar con su ejecución.')
    pause
elseif length(LD.deltae)==1
    LD.deltae(2) = nonExist;
    LD.flags.deltae.dim = 0; %Width==1, no se puede usar pre-lookup
else
    LD.flags.deltae.dim = 1;
end


%% CHECK DELTA OF RUDDER
clear LD.flags.deltar
if isempty(LD.deltar)
    wrn = msgbox('Vector de angulos de ataque vacios: LD.deltar==[] ', 'Aviso','warn');uiwait(wrn);
    disp('Presiona Ctrl+C para detener el programa, o cualquier otra tecla para continuar con su ejecución.')
    pause
elseif length(LD.deltar)==1
    LD.deltar(2) = nonExist;
    LD.flags.deltar.dim = 0; %Width==1, no se puede usar pre-lookup
else
    LD.flags.deltar.dim = 1;
end


%% CHECK RIGHT FLAPERON
clear LD.flags.deltafr
if isempty(LD.deltafr)
    wrn = msgbox('Vector de angulos de ataque vacios: LD.deltafr==[] ', 'Aviso','warn');uiwait(wrn);
    disp('Presiona Ctrl+C para detener el programa, o cualquier otra tecla para continuar con su ejecución.')
    pause
elseif length(LD.deltafr)==1
    LD.deltafr(2) = nonExist;
    LD.flags.deltafr.dim = 0; %Width==1, no se puede usar pre-lookup
else
    LD.flags.deltafr.dim = 1;
end


%% CHECK LEFT FLAPERON
clear LD.flags.deltafl
if isempty(LD.deltafl)
    wrn = msgbox('Vector de angulos de ataque vacios: LD.deltafl==[] ', 'Aviso','warn');uiwait(wrn);
    disp('Presiona Ctrl+C para detener el programa, o cualquier otra tecla para continuar con su ejecución.')
    pause
elseif length(LD.deltafl)==1
    LD.deltafl(2) = nonExist;
    LD.flags.deltafl.dim = 0; %Width==1, no se puede usar pre-lookup
else
    LD.flags.deltafl.dim = 1;
end

