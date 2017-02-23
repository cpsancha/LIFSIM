% rpm=13000*(2*pi)/60;
% Cp=0.0181;
% PWR=4.484; %HP

% rpm=13000*(2*pi)/60;
% Cp=0.0152;
% PWR=3.767; %HP

% rpm=13000*(2*pi)/60;
% Cp=0.0485;
% PWR=11.982; %HP

% rpm=5000/60;
% Cp=0.0272;
% PWR=0.382; %HP

rpm=5000/60;
Cp=0.0272;
PWR=0.382; %HP
J=0.37;
vmph=29.6;
TLbf=3.404;

Mph2Ms   = 0.44704;
Nautical=1852/3600;
vms=vmph*Mph2Ms
vmsNautical=Nautical*vmph

Lbf2N=4.4482216;
TN=TLbf*Lbf2N

D=17*0.0254; %m
rho=1.225;

Hp2Watts = 745.7;

PWRW=PWR*Hp2Watts

% rho=PWRW/(rpm^3*D^5*Cp)
% P=rho*(rpm^3)*(D^5)*Cp
% c=PWRW/P

Vj=J*rpm*D