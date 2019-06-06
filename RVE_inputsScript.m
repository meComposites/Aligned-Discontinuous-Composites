%% Inputs

Nfi=8; Nfj=5;

%Physical inputs
Vf=0.3;
nc=0.5; % number ratio of C over C+G
Lf=3.0;

Dc=0.010; Dg=0.007;

EgAvg=80000; EgCoV=0.07; 
EcAvg=230000; EcCoV=0.05;

XgAvg=3500; XgCoV=0.25;
XcAvg=4500; XcCoV=0.20;

SIAvg=80; SICoV=0.15;
GIAvg=1500; GiicIAvg=0.8;

% Grid of applied strains
nEpsilonRVEGrid=100;
EpsilonRVEMax=0.10;

opts = optimset('TolX',10^(-32),'Display','off');
