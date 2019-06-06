function [xInt,yInt,EInt,deInt] = SPSL_VariableGiicSm_EffectsOfRandomness(Ei,Ti,tmi,Lall,Giicm,Sm,gm)   % Only for single L

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T,tm,Eb,Nx,Ng] = SPInputs_VariableGiicSm(Ei,Ti,tmi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Preliminary calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,type,Gm,y,gm,Sm,Sgm] = SPPreCalc(gm,Sm,Giicm,T,tm,Eb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Matrix of Length of Process Zones %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lpz = SPLpzMatrix(N,type,y,gm,Sm,T,tm,Eb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Loop for all L required %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(Lall)>1
    Nx=1;       % If there are several L, it will not output fields.
end

[x,DX,S,g,nMax,XeRemoteAll,PropertiesAll,XePointsAll,ijActive,DXL,gCrit] = ...
    SPSLLoop(N,type,Gm,y,gm,Sm,Lpz,Lall,T,tm,Eb,Nx,Ng);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Calculates fields if there is only on value of L %%%%%%%%%%%%%%%%%%%%

if length(Lall)==1
    [fieldx,fieldg,fieldS,fieldXB] = ...
        SPPostProcFields(Nx,nMax,x,DX,S,g,DXL);
end

Epsilon=XeRemoteAll(:,1);
Sigma=XeRemoteAll(:,2);

% calculating depsilon
if max(Epsilon)-Epsilon(end)==0
    [~,im]=max(Sigma);
    depsilon=Epsilon(end)-Epsilon(im);
else
    depsilon=Epsilon(end)-max(Epsilon);
end

deltaSigma=3; % (MPa)
Npercents=1.5;  % (%)

% Correction for plateau
[M,indBB]=max(Sigma); % Sigma Max
indMin=min(find(Sigma>max(Sigma)*(1-Npercents/100)));
emin=Epsilon(indMin); % epsilon min in N%
smin=max(Sigma)*(1-Npercents/100);
smin=Sigma(indMin);
indMax=max(find(Sigma>max(Sigma)*(1-Npercents/100)));
emax=max(max(Epsilon(indMin:indMax)),Epsilon(indMax)); % epsilon max in N%
smax=max(Sigma);

% New curve with loading only
if smin==smax
    NewSigma=[Sigma(1:indMin-1) ; smax];
    NewEpsilon=[Epsilon(1:indMin-1) ; emax];
else
    NewSigma=[Sigma(1:indMin-1) ; smin ; smax];
    NewEpsilon=[Epsilon(1:indMin-1) ; emin ; emax];
end

% Interpolation on sigma
NewSigmaGrid=[0:deltaSigma:smax]';

[NewEpsilon,indUnique]=unique(floor(10^5*NewEpsilon)/10^5);
NewSigma=NewSigma(indUnique);

NewEpsilonGrid = interp1(NewSigma,NewEpsilon,NewSigmaGrid);
NewEpsilonGrid(end)=emax;

% % Add unloading depending on toughnessclo
% GicTemp=trapz(NewEpsilonGrid,NewSigmaGrid);
%
%
% dG=Giicm/(Ti)*100-GicTemp;
% depsilon=2*dG/NewSigmaGrid(end);

xInt=NewEpsilonGrid;
yInt=NewSigmaGrid;
EInt=PropertiesAll(2);
deInt=depsilon;


end