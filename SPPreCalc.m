function[N,type,Gm,y,gm,Sm,Sgm] = SPPreCalc(gmIn,SmIn,Giicm,T,tm,Eb)
% Calculates characteristics for the constitutive law

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mortar constitutive law %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add zero point (initial) to mortar constitutive law
gm=[0;gmIn];
Sm=[0;SmIn];

% Calculating number of matrix subdomains, N
N=length(gm);

% Calculating maximum strain, \gamma^[N]
U=(Sm(1:end-1)+Sm(2:end))/2.*(gm(2:end)-gm(1:end-1));           % Eq.8
gmN=gm(end)+2/Sm(end)*(Giicm/tm-sum(U));                        % Eq.8

% Matrix constitutive law, gm=\gamma^[i] vs Sm=\tau^[i], i={0...N+1}
gm=[gm;gmN];
Sm=[Sm;0];
Sgm=[gm*100,Sm];    % strains in percentage


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subdomain parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating tangent modulus, Gm=G^[i]
Gm=(Sm(2:end)-Sm(1:end-1))./(gm(2:end)-gm(1:end-1));

% Defining types of domains
type=zeros(N,1);        % 0 for zero stiffness subdomain
type(Gm>0)=1;           % 1 for positive stiffness subdomain
type(Gm<0)=-1;          % -1 for negative stiffness subdomain

Gm=abs(Gm);

% Characteristic variable for governing differential equation
y=sqrt(2*Gm/(T*tm*Eb));             % Eq. 5, y=\lambda

end