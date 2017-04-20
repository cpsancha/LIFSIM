function stabilityState = ndHeSs2HeSs(ndSs,he,atm,H,varargin)
%ndHeSs2HeSs  transforms a non dimensional stability state in a 
%dimensional one
%FIXME: Polish the function and complete the help.

% Setup options
options = parseOptions(varargin,@setHeroesRigidOptions);

if iscell(he)
    % Dimensions of the output ndHe
    n              = numel(he);
    s              = size(he);
    stabilityState = cell(s);

    % Loop using linear indexing
    for i = 1:n
        stabilityState{i}    = ndHeSs2HeSs_i(ndSs{i},he{i},atm,H,options);
    end

else
    stabilityState    = ndHeSs2HeSs_i(ndSs,he,atm,H,options);
end

end





function Ss = ndHeSs2HeSs_i(ndSs,he,atm,H,varargin)

options   = parseOptions(varargin,@setHeroesRigidOptions);

%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

R        = he.mainRotor.R;

ndTs     = ndSs.ndTs;

Ts       = ndHeTrimState2HeTrimState(ndTs,he,atm,H,options);

ndA      = ndSs.ndA;
ndB      = ndSs.ndB;

ndSD     = ndSs.staDerMatrix.AllElements;
ndSDSum  = ndSs.staDerMatrix.AllElementsFromSum;
ndSDmr   = ndSs.staDerMatrix.mainRotor;
ndSDtr   = ndSs.staDerMatrix.tailRotor;
ndSDf    = ndSs.staDerMatrix.fuselage;
ndSDvf   = ndSs.staDerMatrix.verticalFin;
ndSDlHTP = ndSs.staDerMatrix.leftHTP;
ndSDrHTP = ndSs.staDerMatrix.rightHTP;

ndCD     = ndSs.conDerMatrix.AllElements;
ndCDSum  = ndSs.conDerMatrix.AllElementsFromSum;
ndCDmr   = ndSs.conDerMatrix.mainRotor;
ndCDtr   = ndSs.conDerMatrix.tailRotor;
ndCDf    = ndSs.conDerMatrix.fuselage;
ndCDvf   = ndSs.conDerMatrix.verticalFin;
ndCDlHTP = ndSs.conDerMatrix.leftHTP;
ndCDrHTP = ndSs.conDerMatrix.rightHTP;

Theta   = ndTs.solution.Theta;

sizeTs  = size(Theta);
nFC     = numel(Theta);

%--------------------------------------------------------------------------
% Allocation of variables
%--------------------------------------------------------------------------

A      = zeros([9 9 sizeTs]);
ALO    = zeros([4 4 sizeTs]);
A12    = zeros([4 5 sizeTs]);
A21    = zeros([5 4 sizeTs]);
ALA    = zeros([5 5 sizeTs]);


Atrk   = zeros([8 8 nFC]);
ALOtrk = zeros([4 4 nFC]);
ALAtrk = zeros([4 4 nFC]);

B      = zeros([9 4 sizeTs]);
BLO    = zeros([4 2 sizeTs]);
B12    = zeros([4 2 sizeTs]);
B21    = zeros([5 2 sizeTs]);
BLA    = zeros([5 2 sizeTs]);

SD     = zeros([6 9 sizeTs]);
SDsum  = zeros([6 9 sizeTs]);
SDmr   = zeros([6 9 sizeTs]);
SDtr   = zeros([6 9 sizeTs]);
SDf    = zeros([6 9 sizeTs]);
SDvf   = zeros([6 9 sizeTs]);
SDlHTP = zeros([6 9 sizeTs]);
SDrHTP = zeros([6 9 sizeTs]);

CD     = zeros([6 4 sizeTs]);
CDsum  = zeros([6 4 sizeTs]);
CDmr   = zeros([6 4 sizeTs]);
CDtr   = zeros([6 4 sizeTs]);
CDf    = zeros([6 4 sizeTs]);
CDvf   = zeros([6 4 sizeTs]);
CDlHTP = zeros([6 4 sizeTs]);
CDrHTP = zeros([6 4 sizeTs]);

eigW    = zeros([9 9 sizeTs]);
eigWLO  = zeros([4 4 sizeTs]);
eigWLA  = zeros([5 5 sizeTs]);

si      = zeros([9 sizeTs]);
siLO    = zeros([4 sizeTs]);
siLA    = zeros([5 sizeTs]);

s1 = zeros(sizeTs);
s2 = zeros(sizeTs);
s3 = zeros(sizeTs);
s4 = zeros(sizeTs);
s5 = zeros(sizeTs);
s6 = zeros(sizeTs);
s7 = zeros(sizeTs);
s8 = zeros(sizeTs);
s9 = zeros(sizeTs);

s1LO = zeros(sizeTs);
s2LO = zeros(sizeTs);
s3LO = zeros(sizeTs);
s4LO = zeros(sizeTs);

s1LA = zeros(sizeTs);
s2LA = zeros(sizeTs);
s3LA = zeros(sizeTs);
s4LA = zeros(sizeTs);
s5LA = zeros(sizeTs);

[SDnames,CDnames] = SCDnames;

for s = 1:length(SDnames);
   
        StaDer.(SDnames{s}) = zeros(sizeTs);
      StaDerMr.(SDnames{s}) = zeros(sizeTs);
      StaDerTr.(SDnames{s}) = zeros(sizeTs);
       StaDerF.(SDnames{s}) = zeros(sizeTs);
      StaDerVf.(SDnames{s}) = zeros(sizeTs);
    StaDerLHTP.(SDnames{s}) = zeros(sizeTs);
    StaDerRHTP.(SDnames{s}) = zeros(sizeTs);
    
     StaDerSum.(SDnames{s}) = zeros(sizeTs);
    
end

for c = 1:length(CDnames);
   
        ConDer.(CDnames{c}) = zeros(sizeTs);
      ConDerMr.(CDnames{c}) = zeros(sizeTs);
      ConDerTr.(CDnames{c}) = zeros(sizeTs);
       ConDerF.(CDnames{c}) = zeros(sizeTs);
      ConDerVf.(CDnames{c}) = zeros(sizeTs);
    ConDerLHTP.(CDnames{c}) = zeros(sizeTs);
    ConDerRHTP.(CDnames{c}) = zeros(sizeTs);
    
     ConDerSum.(CDnames{c}) = zeros(sizeTs);
    
end


for i = 1:nFC
    
    
    %----------------------------------------------------------------------
    % Extract nondimensional stability matrix ndA F(1:9)
    %----------------------------------------------------------------------
    
    Omega = Ts.solution.Omega(i);
    
    A(:,:,i)   = ndA2A(ndA(:,:,i),Omega,R);

    ALO(:,:,i) = A((1:4),(1:4),i);
    ALA(:,:,i) = A((5:9),(5:9),i);
    
    A12(:,:,i) = A((1:4),(5:9),i);
    A21(:,:,i) = A((5:9),(1:4),i);
    
    % Stability matrix in linear order for eigenvalues tranking
    Atrk(:,:,i)   = A((1:8),(1:8),i);
    ALOtrk(:,:,i) = A((1:4),(1:4),i);
    ALAtrk(:,:,i) = A((5:8),(5:8),i);

    %----------------------------------------------------------------------
    % Matrix of complete vehicle and components stability derivatives
    %----------------------------------------------------------------------
        
    SD(:,:,i)     = ndStaDer2StaDer(ndSD(:,:,i),Omega,R);
    SDsum(:,:,i)  = ndStaDer2StaDer(ndSDSum(:,:,i),Omega,R);
    SDmr(:,:,i)   = ndStaDer2StaDer(ndSDmr(:,:,i),Omega,R);
    SDtr(:,:,i)   = ndStaDer2StaDer(ndSDtr(:,:,i),Omega,R);
    SDf(:,:,i)    = ndStaDer2StaDer(ndSDf(:,:,i),Omega,R);
    SDvf(:,:,i)   = ndStaDer2StaDer(ndSDvf(:,:,i),Omega,R);
    SDlHTP(:,:,i) = ndStaDer2StaDer(ndSDlHTP(:,:,i),Omega,R);
    SDrHTP(:,:,i) = ndStaDer2StaDer(ndSDrHTP(:,:,i),Omega,R);
    
        
    %----------------------------------------------------------------------
    % Extract nondimensional control matrix ndB F(1:9)
    %----------------------------------------------------------------------
    
    B(:,:,i)   = ndB2B(ndB(:,:,i),Omega,R);
       
    BLO(:,:,i) = B((1:4),(1:2),i);
    BLA(:,:,i) = B((5:9),(3:4),i);
    B12(:,:,i) = B((1:4),(3:4),i);
    B21(:,:,i) = B((5:9),(1:2),i);

    %----------------------------------------------------------------------
    % Matrix of complete vehicle and elementes control derivatives
    %----------------------------------------------------------------------

    CD(:,:,i)     = ndConDer2ConDer(ndCD(:,:,i),Omega,R);
    CDsum(:,:,i)  = ndConDer2ConDer(ndCDSum(:,:,i),Omega,R);
    CDmr(:,:,i)   = ndConDer2ConDer(ndCDmr(:,:,i),Omega,R);
    CDtr(:,:,i)   = ndConDer2ConDer(ndCDtr(:,:,i),Omega,R);
    CDf(:,:,i)    = ndConDer2ConDer(ndCDf(:,:,i),Omega,R);
    CDvf(:,:,i)   = ndConDer2ConDer(ndCDvf(:,:,i),Omega,R);
    CDlHTP(:,:,i) = ndConDer2ConDer(ndCDlHTP(:,:,i),Omega,R);
    CDrHTP(:,:,i) = ndConDer2ConDer(ndCDrHTP(:,:,i),Omega,R);
    
                      
    %----------------------------------------------------------------------
    % Eigen values and eigen vectors calulations
    %----------------------------------------------------------------------
     
    % Solution of the eigenvalues/eigen vector problem for the comple 
    % stability matrix A and the LO submatrix and LA submatrix
    [eigW(:,:,i),dM]     = eig(A(:,:,i));
    [eigWLO(:,:,i),dMLO] = eig(ALO(:,:,i));
    [eigWLA(:,:,i),dMLA] = eig(ALA(:,:,i));
    
     % eigenvalues indexing and structuring 
     si(:,i)   = diag(dM);
     siLO(:,i) = diag(dMLO);
     siLA(:,i) = diag(dMLA);
    
     s1(i) = si(1,i);
     s2(i) = si(2,i);
     s3(i) = si(3,i);
     s4(i) = si(4,i);
     s5(i) = si(5,i);
     s6(i) = si(6,i);
     s7(i) = si(7,i);
     s8(i) = si(8,i);
     s9(i) = si(9,i);
    
     s1LO(i) = siLO(1,i);
     s2LO(i) = siLO(2,i);
     s3LO(i) = siLO(3,i);
     s4LO(i) = siLO(4,i);
     
     s1LA(i) = siLA(1,i);
     s2LA(i) = siLA(2,i);
     s3LA(i) = siLA(3,i);
     s4LA(i) = siLA(4,i);
     s5LA(i) = siLA(5,i);
     
    %---------------------------------------------------------------------
          
    %---------------------------------------------------------------------- 
    % Generation of substruct ndStaDer with non dimensional stability
    % derivatives
    %----------------------------------------------------------------------
    
    for s = 1:length(SDnames);
   
        matrixSD = (SD(:,:,i))';
      matrixSDmr = (SDmr(:,:,i))';
      matrixSDtr = (SDtr(:,:,i))';
       matrixSDf = (SDf(:,:,i))';
      matrixSDvf = (SDvf(:,:,i))';
    matrixSDlHTP = (SDlHTP(:,:,i))';
    matrixSDrHTP = (SDrHTP(:,:,i))';
    matrixSDsum  = (SDsum(:,:,i))';    
    
    
        StaDer.(SDnames{s})(i) = matrixSD(s);
      StaDerMr.(SDnames{s})(i) = matrixSDmr(s);
      StaDerTr.(SDnames{s})(i) = matrixSDtr(s);
       StaDerF.(SDnames{s})(i) = matrixSDf(s);
      StaDerVf.(SDnames{s})(i) = matrixSDvf(s);
    StaDerLHTP.(SDnames{s})(i) = matrixSDlHTP(s);
    StaDerRHTP.(SDnames{s})(i) = matrixSDrHTP(s);
     StaDerSum.(SDnames{s})(i) = matrixSDsum(s);
       
    end 
    %----------------------------------------------------------------------
     
    %---------------------------------------------------------------------- 
    %Generation of substruct ndStaDer with non dimensional control
    %derivatives
    %---------------------------------------------------------------------- 
    
    for c = 1:length(CDnames);
   
        matrixCD = (CD(:,:,i))';
      matrixCDmr = (CDmr(:,:,i))';
      matrixCDtr = (CDtr(:,:,i))';
       matrixCDf = (CDf(:,:,i))';
      matrixCDvf = (CDvf(:,:,i))';
    matrixCDlHTP = (CDlHTP(:,:,i))';
    matrixCDrHTP = (CDrHTP(:,:,i))';   
     matrixCDsum = (CDsum(:,:,i))';    
    
      ConDer.(CDnames{c})(i) = matrixCD(c);
    ConDerMr.(CDnames{c})(i) = matrixCDmr(c);
    ConDerTr.(CDnames{c})(i) = matrixCDtr(c);
     ConDerF.(CDnames{c})(i) = matrixCDf(c);
    ConDerVf.(CDnames{c})(i) = matrixCDvf(c);
  ConDerLHTP.(CDnames{c})(i) = matrixCDlHTP(c);
  ConDerRHTP.(CDnames{c})(i) = matrixCDrHTP(c);
   ConDerSum.(CDnames{c})(i) = matrixCDsum(c);
    
    end
    
%----------------------------------------------------------------------
    
end

%==========================================================================
% Tracked eigen values and eigen vectors calulations
%==========================================================================

%--------------------------------------------------------------------------
% Solution of the tracked eigenvalues and eigenvector problem
%--------------------------------------------------------------------------
 [Wtrack,siTrack]     = eigenshuffle(Atrk);
 [WLOtrack,siLOTrack] = eigenshuffle(ALOtrk);
 [WLAtrack,siLATrack] = eigenshuffle(ALAtrk);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Trakeg eigenvalues and eigenvectors allocation
%--------------------------------------------------------------------------
 eigWtr    = zeros([8 8 sizeTs]);
 eigWLOtr  = zeros([4 4 sizeTs]);
 eigWLAtr  = zeros([4 4 sizeTs]);
 
 siTr      = zeros([8 sizeTs]);
 siLOtr    = zeros([4 sizeTs]);
 siLAtr    = zeros([4 sizeTs]);
 
 s1Tr      = zeros(sizeTs);
 s2Tr      = zeros(sizeTs);
 s3Tr      = zeros(sizeTs);
 s4Tr      = zeros(sizeTs);
 s5Tr      = zeros(sizeTs);
 s6Tr      = zeros(sizeTs);
 s7Tr      = zeros(sizeTs);
 s8Tr      = zeros(sizeTs);
 s9Tr      = zeros(sizeTs);

 s1LOtr    = zeros(sizeTs);
 s2LOtr    = zeros(sizeTs);
 s3LOtr    = zeros(sizeTs);
 s4LOtr    = zeros(sizeTs);

 s1LAtr    = zeros(sizeTs);
 s2LAtr    = zeros(sizeTs);
 s3LAtr    = zeros(sizeTs);
 s4LAtr    = zeros(sizeTs);
 s5LAtr    = zeros(sizeTs);
 %-------------------------------------------------------------------------

 %=========================================================================
 % Main loop for tracked eigen values structuring
 %=========================================================================
 
for i = 1:nFC;
    
 siTr(:,i)   = siTrack(:,i);
 siLOtr(:,i) = siLOTrack(:,i);
 siLAtr(:,i) = siLATrack(:,i);
    
 s1Tr(i)     = siTrack(1,i);
 s2Tr(i)     = siTrack(2,i);
 s3Tr(i)     = siTrack(3,i);
 s4Tr(i)     = siTrack(4,i);
 s5Tr(i)     = siTrack(5,i);
 s6Tr(i)     = siTrack(6,i);
 s7Tr(i)     = siTrack(7,i);
 s8Tr(i)     = siTrack(8,i);
 s9Tr(i)     = 0;

 s1LOtr(i)   = siLOTrack(1,i);
 s2LOtr(i)   = siLOTrack(2,i);
 s3LOtr(i)   = siLOTrack(3,i);
 s4LOtr(i)   = siLOTrack(4,i);

 s1LAtr(i)   = siLATrack(1,i);
 s2LAtr(i)   = siLATrack(2,i);
 s3LAtr(i)   = siLATrack(3,i);
 s4LAtr(i)   = siLATrack(4,i);
 s5LAtr(i)   = 0;
    
 eigWtr(:,:,i)   = Wtrack(:,:,i);
 eigWLOtr(:,:,i) = WLOtrack(:,:,i);
 eigWLAtr(:,:,i) = WLAtrack(:,:,i);
 
end

%--------------------------------------------------------------------------
% Calculation of characteristic values associated to tracked eigen values
% FIXME: ALVARO, We must make something more consistent depending
% on the nature of the real and imag part of the eigenvalues
%--------------------------------------------------------------------------

 modsiTr    = abs(siTr);
 omegasiTr  = abs(imag(siTr));
 zetasiTr   = -real(siTr)./abs(imag(siTr));
 omega0siTr = omegasiTr./sqrt(1-zetasiTr.^2);
 invTausiTr = 1./real(siTr);
 t1_2siTr   = log(1/2)*(real(siTr)).^(-1);
 t2siTr     = log(2)*(real(siTr)).^(-1);
 
 modsiLOtr    = abs(siLOtr);
 omegasiLOtr  = abs(imag(siLOtr));
 zetasiLOtr   = -real(siLOtr)./abs(imag(siLOtr));
 omega0siLOtr = omegasiLOtr./sqrt(1-zetasiLOtr.^2);
 invTausiLOtr = 1./real(siLOtr);
 t1_2siLOtr   = log(1/2)*(real(siLOtr)).^(-1);
 t2siLOtr     = log(2)*(real(siLOtr)).^(-1);
 
 modsiLAtr    = abs(siLAtr);
 omegasiLAtr  = abs(imag(siLAtr));
 zetasiLAtr   = -real(siLAtr)./abs(imag(siLAtr));
 omega0siLAtr = omegasiLAtr./sqrt(1-zetasiLAtr.^2);
 invTausiLAtr = 1./real(siLAtr);
 t1_2siLAtr   = log(1/2)*(real(siLAtr)).^(-1);
 t2siLAtr     = log(2)*(real(siLAtr)).^(-1);

%==========================================================================

%==========================================================================
% Structuring the results
%==========================================================================
      
    
    
staDerMatrix = struct('AllElements',SD,'AllElementsFromSum',SDsum,...
                 'mainRotor',SDmr,'tailRotor',SDtr,...
                 'fuselage',SDf,'verticalFin',SDvf,...
                 'leftHTP',SDlHTP,'rightHTP',SDrHTP);
             
conDerMatrix = struct('AllElements',CD,'AllElementsFromSum',CDsum,...
                 'mainRotor',CDmr,'tailRotor',CDtr,...
                 'fuselage',CDf,'verticalFin',CDvf,...
                 'leftHTP',CDlHTP,'rightHTP',CDrHTP);             

staDer = struct('AllElements',StaDer,'AllElementsFromSum',StaDerSum,...
                 'mainRotor',StaDerMr,'tailRotor',StaDerTr,...
                 'fuselage',StaDerF,'verticalFin',StaDerVf,...
                 'leftHTP',StaDerLHTP,'rightHTP',StaDerRHTP);
             
conDer = struct('AllElements',ConDer,'AllElementsFromSum',ConDerSum,...
                 'mainRotor',ConDerMr,'tailRotor',ConDerTr,...
                 'fuselage',ConDerF,'verticalFin',ConDerVf,...
                 'leftHTP',ConDerLHTP,'rightHTP',ConDerRHTP);
             
eigenVal = struct('s1',s1,'s2',s2,'s3',s3,'s4',s4,'s5',s5,...
                   's6',s6,'s7',s7,'s8',s8,'s9',s9);
               
eigenValLO = struct('s1LO',s1LO,'s2LO',s2LO,'s3LO',s3LO,'s4LO',s4LO);  
 
eigenValLA = struct('s1LA',s1LA,'s2LA',s2LA,'s3LA',s3LA,...
              's4LA',s4LA,'s5LA',s5LA);
          
           
eigenValTr = struct('s1',s1Tr,'s2',s2Tr,'s3',s3Tr,'s4',s4Tr,'s5',s5Tr,...
                   's6',s6Tr,'s7',s7Tr,'s8',s8Tr,'s9',s9Tr);
               
eigenValLOtr = struct('s1LO',s1LOtr,'s2LO',s2LOtr,...
                       's3LO',s3LOtr,'s4LO',s4LOtr);  
 
eigenValLAtr = struct('s1LA',s1LAtr,'s2LA',s2LAtr,...
                       's3LA',s3LAtr,'s4LA',s4LAtr,'s5LA',s5LAtr);
                   

charValsTr = struct('mod',modsiTr,'omega',omegasiTr,'zeta',zetasiTr,...
                  'omegaN',omega0siTr,'invTau',invTausiTr,...
                  't1_2',t1_2siTr,'t2',t2siTr);
              
charValsLOTr = struct('mod',modsiLOtr,'omega',omegasiLOtr,'zeta',zetasiLOtr,...
                  'omegaN',omega0siLOtr,'invTau',invTausiLOtr,...
                  't1_2',t1_2siLOtr,'t2',t2siLOtr);
              
charValsLATr = struct('mod',modsiLAtr,'omega',omegasiLAtr,'zeta',zetasiLAtr,...
                  'omegaN',omega0siLAtr,'invTau',invTausiLAtr,...
                  't1_2',t1_2siLAtr,'t2',t2siLAtr);
                   

stabilityDerivatives = struct(...
             'A',A, ...   
             'ALO',ALO, ...
             'A12',A12, ...
             'A21',A21, ...
             'ALA',ALA, ...
             'staDerMatrix',staDerMatrix, ...
             'staDer',staDer  ...
);

controlDerivatives   = struct(...
             'B',B, ...
             'BLO',BLO, ...
             'B12',B12, ...
             'B21',B21, ...
             'BLA',BLA, ...
             'conDerMatrix',conDerMatrix, ...
             'conDer',conDer  ...
);
eigenSolution        = struct(...
             'eigW', eigW, ...
             'eigWLO',eigWLO, ...
             'eigWLA',eigWLA, ...
             'si',si, ...
             'siLO',siLO, ...
             'siLA',siLA, ...
             'eigenVal',eigenVal, ...
             'eigenValLO',eigenValLO, ...
             'eigenValLA',eigenValLA,...
             'eigWtr', eigWtr, ...
             'eigWLOtr',eigWLOtr, ...
             'eigWLAtr',eigWLAtr, ...
             'siTr',siTr, ...
             'siLOtr',siLOtr, ...
             'siLAtr',siLAtr, ...
             'eigenValTr',eigenValTr, ...
             'eigenValLOtr',eigenValLOtr, ...
             'eigenValLAtr',eigenValLAtr,...
             'charValTr',charValsTr,...
             'charValLOTr',charValsLOTr,...
             'charValLATr',charValsLATr ...
);

Ss                   = struct(...
             'stabilityDerivatives',stabilityDerivatives,...
             'controlDerivatives',controlDerivatives,...
             'eigenSolution',eigenSolution ...
);

% % % % % % % % % % % Ss  = struct('Ts',Ts, ...
% % % % % % % % % % %              'A',A, ...   
% % % % % % % % % % %              'B',B, ...
% % % % % % % % % % %              'ALO',ALO, ...
% % % % % % % % % % %              'A12',A12, ...
% % % % % % % % % % %              'A21',A21, ...
% % % % % % % % % % %              'ALA',ALA, ...
% % % % % % % % % % %              'BLO',BLO, ...
% % % % % % % % % % %              'B12',B12, ...
% % % % % % % % % % %              'B21',B21, ...
% % % % % % % % % % %              'BLA',BLA, ...
% % % % % % % % % % %              'staDerMatrix',staDerMatrix, ...
% % % % % % % % % % %              'conDerMatrix',conDerMatrix, ...
% % % % % % % % % % %              'staDer',staDer, ...
% % % % % % % % % % %              'conDer',conDer, ...
% % % % % % % % % % %              'eigW', eigW, ...
% % % % % % % % % % %              'eigWLO',eigWLO, ...
% % % % % % % % % % %              'eigWLA',eigWLA, ...
% % % % % % % % % % %              'si',si, ...
% % % % % % % % % % %              'siLO',siLO, ...
% % % % % % % % % % %              'siLA',siLA, ...
% % % % % % % % % % %              'eigenVal',eigenVal, ...
% % % % % % % % % % %              'eigenValLO',eigenValLO, ...
% % % % % % % % % % %              'eigenValLA',eigenValLA,...
% % % % % % % % % % %              'eigWtr', eigWtr, ...
% % % % % % % % % % %              'eigWLOtr',eigWLOtr, ...
% % % % % % % % % % %              'eigWLAtr',eigWLAtr, ...
% % % % % % % % % % %              'siTr',siTr, ...
% % % % % % % % % % %              'siLOtr',siLOtr, ...
% % % % % % % % % % %              'siLAtr',siLAtr, ...
% % % % % % % % % % %              'eigenValTr',eigenValTr, ...
% % % % % % % % % % %              'eigenValLOtr',eigenValLOtr, ...
% % % % % % % % % % %              'eigenValLAtr',eigenValLAtr,...
% % % % % % % % % % %              'charValTr',charValsTr,...
% % % % % % % % % % %              'charValLOTr',charValsLOTr,...
% % % % % % % % % % %              'charValLATr',charValsLATr ...
% % % % % % % % % % %              );

end

%==========================================================================
% TO BE INCLUDED IN THE HELP
%==========================================================================

%              u  |  w  |  oy  |  Th   ||  v  |  ox  |  Ph  |  ox  |  Ps  | 
%(du/dt) :F1 |    |     |      |       ||     |      |      |      |      | 
%(dw/dt) :F2 |
%(doy/dt):F3 |
%(dTh/dt):F4 |
%
%(dv/dt) :F5 |
%(dox/dt):F6 |
%(dPh/dt):F7 |
%(doz/dt):F8 |
%(dPs/dt):F9 |
%
%
%On the other hand, the control matrix ndB
%
%               t0  |  t1S  |  t1C  |  t0tr  |  
%(du/dt) :F1 |              
%(dw/dt) :F2 |
%(doy/dt):F3 |
%(dTh/dt):F4 |
%
%(dv/dt) :F5 |
%(dox/dt):F6 |
%(dPh/dt):F7 |
%(doz/dt):F8 |
%(dPs/dt):F9 |


