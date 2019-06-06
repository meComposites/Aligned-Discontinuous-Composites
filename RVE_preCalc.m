function [SumEv,DeltaEv,SumEh,DeltaEh,Eneighbours,Eav,Ebv,Eah,Ebh,...
    statusv0,statush0,xDiscf,nDiscf,...
    EpsilonRVE] = RVE_preCalc(Ef,Nfi,Nfj,EpsilonRVEMax,nEpsilonRVEGrid,xEndf)
% Preliminary calculations for the RVE 
%(constants throughout the entire analysis)
%and initiation of variables which are changed

%% Calculation & mapping of moduli
[Eav,Ebv,Eah,Ebh] = aux_mapping_f_to_I(Ef);

[SumEv,DeltaEv] = aux_SumDelta(Eav,Ebv);
[SumEh,DeltaEh] = aux_SumDelta(Eah,Ebh);

% Map modulus: neighbours from fibres
Eneighbours = inf(Nfi,Nfj,4);
%1, above:
Eneighbours(2:end,:,1)=Ef(1:end-1,:);
%2, on the right:
Eneighbours(:,1:end-1,2)=Ef(:,2:end);
%3, below:
Eneighbours(1:end-1,:,3)=Ef(2:end,:);
%4, on the left
Eneighbours(:,2:end,4)=Ef(:,1:end-1);

%% Initiating interactions 
%Sets the status of all ineractions to "1=needs updating"
statusv0=ones(Nfi-1,Nfj);
statush0=ones(Nfi,Nfj-1);

%% Fibre discontinuities and number of breaks
% 3D Matrix of discontinuities in fibres 
% (cross section * 5 discontinuities maximum) 
xDiscf=nan(Nfi,Nfj,5);
xDiscf(:,:,1)=xEndf;
nDiscf=sum(~isnan(xDiscf),3);
%% Creates RVE strain vector
EpsilonRVE=linspace(0,EpsilonRVEMax,nEpsilonRVEGrid)';

end

