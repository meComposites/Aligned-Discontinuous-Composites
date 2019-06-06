function [SigmaFavg,...
        kfCrit,DeNextBreak,iFBreak,jFBreak] = ...
        Fibres_SS(Nfi,Nfj,Df,Tf,Ef,Eneighbours,...
        Sigmav,Sigmah,XfNow,nEpsilonNow)
% Calculates the SS-curve of the fibres in the RVE,
% and returns the failed fibres and critical neighbours

%% Calculates factors to go from SigmaI to SigmaF (Eq 18 & 19 in paper)

% Normalised area (Tf*Df) for each fibre, for the 4 interactions
TtimesDf=Tf.*Df;

% Normalised area (Tf*Df) for the neighbours of each fibre
TtimesDneighbours=zeros(Nfi,Nfj,4);
%1, above (neighbour 1 has interaction 3 with this fibre)
TtimesDneighbours(2:end,:,1)=Tf(1:end-1,:,3).*Df(1:end-1,:);
%2, on the right (neighbour 2 has interaction 4 with this fibre)
TtimesDneighbours(:,1:end-1,2)=Tf(:,2:end,4).*Df(:,2:end);
%3, below (neighbour 3 has interaction 1 with this fibre)
TtimesDneighbours(1:end-1,:,3)=Tf(2:end,:,1).*Df(2:end,:);
%4, on the left (neighbour 4 has interaction 1 with this fibre)
TtimesDneighbours(:,2:end,4)=Tf(:,1:end-1,2).*Df(:,1:end-1);

% Factors to go from Interaction stresses to Fibre stresses
% (need to avoid NaN created by 0/0 when the neighbours do not exist)
factorSigmaItoFavg=Ef.*(TtimesDf.^2+TtimesDneighbours.^2)./...
    (Ef.*TtimesDf.^2+Eneighbours.*TtimesDneighbours.^2);
factorSigmaItoFavg(isnan(factorSigmaItoFavg))=0;

factorSigmaItoFpeak=(TtimesDf.^2+TtimesDneighbours.^2)./(TtimesDf.^2);
factorSigmaItoFpeak(isnan(factorSigmaItoFpeak))=0;

%% Calculates average fibre stresses from interaction stresses

SigmaIwithNeighbours = zeros(Nfi,Nfj,4,nEpsilonNow);

%1, above:
SigmaIwithNeighbours(2:end,:,1,:) = permute(Sigmav,[1 2 4 3]);
%2, on the right:
SigmaIwithNeighbours(:,1:end-1,2,:) = permute(Sigmah,[1 2 4 3]);
%3, below:
SigmaIwithNeighbours(1:end-1,:,3,:) = permute(Sigmav,[1 2 4 3]);
%4, on the left
SigmaIwithNeighbours(:,2:end,4,:) = permute(Sigmah,[1 2 4 3]);

SigmaFavg=permute(sum(SigmaIwithNeighbours.*factorSigmaItoFavg.*(Tf./Df),3),[1 2 4 3]);

%% Calculates peak fibre stresses and critical neighbour

DeltaSigmaFpeakfromNeighbours=SigmaIwithNeighbours.*...
    (factorSigmaItoFpeak-factorSigmaItoFavg).*(Tf./Df);

% Calculates maximum additional stress peak and corresponding neighbour 
[DeltaSigmaFMax,kfCrit]=max(permute(DeltaSigmaFpeakfromNeighbours,[1 2 4 3]),[],4);

% Checks for failure (if maximum fibre stress is above fibre strength)
fBreakNow=(SigmaFavg+DeltaSigmaFMax)>=XfNow;
DeNextBreak=find(sum(sum(fBreakNow,1),2),1);
[iFBreak,jFBreak]=find(fBreakNow(:,:,DeNextBreak));

end

