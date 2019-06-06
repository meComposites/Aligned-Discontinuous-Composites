function [nc,xEndf,Typef,Df,Ef,Xf,...
    tv,th,Gv,Gh,Sv,Sh,Giicv,Giich]...
    = RVE_generate...
    (Nfi,Nfj,Vf,Lf,ncTarget,...
    Dg,Dc,EgAvg,EgCoV,EcAvg,EcCoV,...
    XgAvg,XgCoV,XcAvg,XcCoV,...
    SIAvg,SICoV,GIAvg,GiicIAvg)

%Calculates the microstructure of the RVE 
%(for the fibres, and for vertical and horizontal interfaces)

%% Fibre variables
%Matrix of fibre-types
[Typef,nc]=defineTypef(Nfi,Nfj,ncTarget);

%Matrix of fibre-ends
xEndf=Lf*rand(Nfi,Nfj);

%Matrix of fibre-diameters
Df=(1-Typef)*Dg+Typef*Dc;

%Matrix of fibre-moduli
Ef=defineEf(Typef,EgAvg,EgCoV,EcAvg,EcCoV);

%Matrix of fibre-strengths
Xf=defineXf(Nfi,Nfj,Typef,XgAvg,XgCoV,XcAvg,XcCoV);

%% Interaction variables
%Matrices of interaction shear strength
[Sv,Sh]=defineSI(Nfi,Nfj,SIAvg,SICoV);

%Matrices of interaction fracture toughness and shear moduli
Gv=GIAvg*ones(Nfi-1,Nfj);
Gh=GIAvg*ones(Nfi,Nfj-1);
Giicv=GiicIAvg*ones(Nfi-1,Nfj);
Giich=GiicIAvg*ones(Nfi,Nfj-1);

%Matrices of interaction thickness
[tv,th]=definetI(Nfi,Nfj,Vf,nc,Dg,Dc);

%Matrix of fibre thicknesses for interactions with k=1:4 neighbours
Tf = defineTfI(Nfi,Nfj,Df,tv,th);

end

%% Defines matrix of fibre-types T
function [Typef,nc] = defineTypef(Nfi,Nfj,ncTarget)

%Nfi = total number of fibres in the y-direction
%Nfj = total number of fibres in the x-direction
%nc = number fraction of C-fibres

%Typef = matrix of fibre-types. Typef=0 for G, Typef=1 for C.

Typef=rand(Nfi,Nfj);
Typef(Typef<(1-ncTarget))=0;  % number-fraction of 1-nc fibres are G
Typef(Typef>0)=1;       % other fibres are C

nc=sum(Typef(:))/(Nfi*Nfj); %update number-fraction of C

end

%% Defines matrix of fibre-types moduli E
function Ef = defineEf(Typef,EgAvg,EgCoV,EcAvg,EcCoV)

%Typef = matrix of fibre-types. Typef=0 for G, Typef=1 for C.
%Requires Expected (Avg) and CoV=Std/Avg values of normal distribution

%randn = from N(0,1). N(Avg,Std)=Avg+N(0,1)*Std=Avg*(1+N(0,1)*CoV)

Ef=(1-Typef)*EgAvg.*(1+randn(size(Typef)).*EgCoV) +...  %if G
    Typef*EcAvg.*(1+randn(size(Typef)).*EcCoV);         %if C

end

%% Defines matrix of fibre-types strengths X
function Xf = defineXf(Nfi,Nfj,Typef,XgAvg,XgCoV,XcAvg,XcCoV)

%Typef = matrix of fibre-types. Typef=0 for G, Typef=1 for C.
%Requires Expected (Avg) and CoV=Std/Avg values of Weibull distribution.
%At the moment, this is NOT scaled to length or stress field!

%Weibull parameters
XgWm=fzero(@(mg) sqrt(gamma(1+2/mg)/(gamma(1+1/mg)^2)-1)-XgCoV, 1.2/XgCoV);
XgWs=XgAvg/(gamma(1+1/XgWm));

XcWm=fzero(@(mc) sqrt(gamma(1+2/mc)/(gamma(1+1/mc)^2)-1)-XcCoV, 1.2/XcCoV);
XcWs=XcAvg/(gamma(1+1/XcWm));

Xf=zeros(Nfi,Nfj,5);

Xf(:,:,1:4)=(1-Typef).*wblrnd(XgWs,XgWm,Nfi,Nfj,4)+...    %if G
    Typef.*wblrnd(XcWs,XcWm,Nfi,Nfj,4);                   %if C

Xf(:,:,1:4)=sort(Xf(:,:,1:4),3);

end

%% Defines matrix of shear lag strength SI
function [Sv,Sh] = defineSI(Nfi,Nfj,SmAvg,SmCoV)

%SmV = matrix of shear strength in vertical (along i) direction.
%SmH = matrix of shear strength in horizontal (along j) direction.
%Requires Expected (Avg) and CoV=Std/Avg values of Weibull distribution

%Weibull parameters
SmWm=fzero(@(m) sqrt(gamma(1+2/m)/(gamma(1+1/m)^2)-1)-SmCoV, 1.2/SmCoV);
SmWs=SmAvg/(gamma(1+1/SmWm));

Sv=wblrnd(SmWs,SmWm,Nfi-1,Nfj);
Sh=wblrnd(SmWs,SmWm,Nfi,Nfj-1);

end

%% Defines matrix of shear lag thicknesses tI
function [tv,th] = definetI(Nfi,Nfj,Vf,nc,Dg,Dc)

NGridtI=1000;

%Calculating average interaction thickness:
DfAvg=(1-nc)*Dg+nc*Dc;
Df2Avg=(1-nc)*Dg^2+nc*Dc^2;

tIAvg = 1024/(405*pi) * DfAvg *...
    ( sqrt(1+ (405*pi)/1024 * Df2Avg/DfAvg^2 * (1-Vf)/Vf) - 1);

%Parameter for distribution (areal density of events)
p=1024/(81*pi^2)/tIAvg^2;

%Generate random numbers
Ftv=rand(Nfi-1,Nfj);
Fth=rand(Nfi,Nfj-1);

%Convert random numbers into random values of tI
%Finds maximum tI required. 
%Initial searching point is 10*tIAvg, and assigns the maximum 
%for a 1.1*larger thickness than the supposed maximum 
%(to avoid extrapolation)
tIMax=fzero(@(tI)...
    1 - (1-pi*p*tI^2/6)*erfc(sqrt(pi*p)*tI/2) ...
    - tI*sqrt(p)*exp(-pi*p*tI^2/4) ...
    -max([Ftv(:);Fth(:)]),10*tIAvg)*1.1;

%Creating grid with CDF:
%Creates grid of tI values
tIGrid=linspace(0,tIMax,NGridtI);
%Calculates CDF at the grid values
CDFtIGrid = 1 - (1-pi*p*tIGrid.^2/6).*erfc(sqrt(pi*p)*tIGrid/2) ...
    - sqrt(p)*tIGrid.*exp(-pi*p*tIGrid.^2/4);

%Calculating tI matrices
tv=interp1(CDFtIGrid,tIGrid,Ftv);
th=interp1(CDFtIGrid,tIGrid,Fth);

end

%% Defines matrix of fibre thicknesses for the interaction with neighbours, TfI
function TfI = defineTfI(Nfi,Nfj,Df,tv,th)

%Calculating distance from fibre(i,j) to its k=1:4 neighbours (3D matrices):
distfNeighbours=inf(Nfi,Nfj,4);
%1, above:
distfNeighbours(2:end,:,1)=Df(1:end-1,:)/2+tv;
%2, on the right:
distfNeighbours(:,1:end-1,2)=Df(:,2:end)/2+th;
%3, below:
distfNeighbours(1:end-1,:,3)=Df(2:end,:)/2+tv;
%4, on the left
distfNeighbours(:,2:end,4)=Df(:,1:end-1)/2+th;

%Define lambda (load sharing coefficient) for each fibre
y=1./sum(1./distfNeighbours.^2,3);

%Calculate fibre thicknesses for the interaction with each neighbour
TfI=repmat(y.*Df,1,1,4).*1./distfNeighbours.^2;

end