function [SigmaIAtRVEGrid] = combineSeriesI(...
    nSLI,lSLI,EpsilonSLI,SigmaSLI,ESLI,deltaeSLI,...
    lBra,lBrb,Ea,Eb,Lf,SumTI,Ta,Tb,...
    EpsilonRVEGrid,DeltaEpsilonRVE)
% Combines all segments of an interaction in series, and interpolates to
% RVE grid of Epsilon

%% Interpolating strains to given stresses
% Step in stress
DSigma=3; % (MPa)

% Strength of the interation corresponds to the first (shortest) SL segment:
[SigmaMaxSL,iStrengthSL]=max(SigmaSLI);

% Lowest secant modulus of all cases
%

%Creating interpolation vector
nSigma=floor(SigmaMaxSL(1)/DSigma);
SigmaIGrid=(0:DSigma:nSigma*DSigma)';

% Initiating matrix of interpolated displacements (strains*segment length):
uSLIAtGrid=nan(nSigma+1,nSLI);

for k=1:nSLI
    uSLIAtGrid(:,k)=lSLI(k)*interp1q(...
        SigmaSLI(1:iStrengthSL(k),k),EpsilonSLI(1:iStrengthSL(k),k),SigmaIGrid);
end

%% Combining in series - up to SigmaMax

% Strains (%) during loading stage
EpsilonIAtGrid=1/Lf*(sum(uSLIAtGrid,2)+...
    SumTI/Tb*SigmaIGrid/Eb*lBra+...
    SumTI/Ta*SigmaIGrid/Ea*lBrb);

%% Combining in series - after SigmaMax

EpsilonIAtGrid(end+1)=EpsilonIAtGrid(end)+max([10^-6 ...
    (deltaeSLI(1)*lSLI(1)/Lf...         %failure of weakest SL 
    -SigmaMaxSL(1)/Lf*...            %minus contibution from elastic unloading:
    (sum(lSLI(2:end)./ESLI(2:end))...   %...of other SL
    +lBra/(2*Eb)+lBrb/(2*Ea)))]);         %...of Br

SigmaIGrid(end+1)=0;

%% Interpolating SigmaI to Grid of RVE Epsilon

%Only interpolates (not extrapolates) the existing part of the curve, up to
%the last calculated value of EpsilonIAtGrid
nPointsForRVEGrid=min(max(0,floor(...
    (EpsilonIAtGrid(end)-EpsilonRVEGrid(1))/DeltaEpsilonRVE)+1),...
    length(EpsilonRVEGrid));

SigmaIAtRVEGrid(1:nPointsForRVEGrid,1) = ...
    interp1q(EpsilonIAtGrid,SigmaIGrid,...
    EpsilonRVEGrid(1:nPointsForRVEGrid,1));

SigmaIAtRVEGrid(nPointsForRVEGrid+1:length(EpsilonRVEGrid),1)=0;

