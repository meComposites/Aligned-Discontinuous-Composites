function [Eini,Lchar,sigma_inf,epsilon_inf] = Eini_2Stiffness(DE,DT,EE,ET,tm,L,G)

% Modified 13 August 2017 SP: 
%   Strains in absolute values (not %)
%   Modulus in MPa (not GPa)
%   Removed commented code

% Definition of E and T
EA=(EE-DE)/2;
EB=(EE+DE)/2;

TA=(ET-DT)/2;
TB=(ET+DT)/2;

% definition of other parameters
K=2*(DE*ET)/(EE*ET+DE*DT);
y=sqrt(8*abs(G(1))/tm*(EE*ET+DE*DT)/((EE^2-DE^2)*(ET^2-DT^2)));
%y2=sqrt(8*abs(G(1))/tm*(EE*ET)/((EE^2)*(ET^2)));
R=(EA*TA)/(EB*TB);

% Characteristic length
Lchar=1/(2*y)*log((R+exp(2*y*L))/(R+exp(-2*y*L)));
% if Lchar == inf     %if 2*y*L too big, Lchar converges towards L
% Lchar=L;
% end

% Stress and strain 
nu=4*ET/(ET^2-DT^2)*1/y;
sigma_inf=nu*sinh(y*Lchar)/(ET/TA+K); % * S0
%sigma_inf2=ET/(y*TB)*sinh(y*Lchar)/(ET+TA*K); % * S0
epsilon_inf=tm/(4*L)*1/G*(cosh(y*Lchar)*(EE+DE)*(ET+DT) + cosh(y*(2*L-Lchar))*(EE-DE)*(ET-DT))/(EE*ET+DE*DT)+2*sigma_inf*ET/(EE*ET+DE*DT); % *S0
if sigma_inf == inf     %if 2*y*L too big, cosh ~ sinh
sigma_inf=nu/(ET/TA+K); % * S0
epsilon_inf=tm/(4*L)*1/G*((EE+DE)*(ET+DT) + (EE-DE)*(ET-DT))/(EE*ET+DE*DT)+2*sigma_inf*ET/(EE*ET+DE*DT); % *S0
end

% Initial stiffness
Eini=sigma_inf/epsilon_inf; % in MPA
Eini=Eini*(ET)/(ET+tm); % account for matrix thickness

Eav=EE/2;
num=1+(DE*DT)/(EE*ET);

E=Eav*num/(1+1/(y*L*tanh(y*Lchar))*(cosh(2*y*L)/2-tanh(y*Lchar)*sinh(2*y*L)/2+1/2/R));
E=E*(ET)/(ET+tm);


end
