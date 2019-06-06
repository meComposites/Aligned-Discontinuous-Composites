function [SigmaRVE,eMax,XRVEMax,iCluster,jCluster]=...
    RVE_SS(~,xEndf,~,Df,Ef,Xf,...
    tv,th,Gv,Gh,Sv,Sh,Giicv,Giich,...
    Lf,nEpsilonRVE,EpsilonRVEMax,~,Vf,nfInCluster,Nfi,Nfj)


%% Pre-calculation of parameters constant throughout the entire analysis
[SumEv,DeltaEv,SumEh,DeltaEh,Eneighbours,Eav,Ebv,Eah,Ebh,...
    statusv,statush,xDiscf,~,...
    EpsilonRVE] = RVE_preCalc(Ef,Nfi,Nfj,EpsilonRVEMax,nEpsilonRVE,xEndf);

XfNow=Xf(:,:,1);
ePrevBreak=0;
SigmaFavg=zeros(Nfi,Nfj,nEpsilonRVE);
SigmavRVE=zeros(Nfi-1,Nfj,nEpsilonRVE);
SigmahRVE=zeros(Nfi,Nfj-1,nEpsilonRVE);

%% Runs loop in EpsilonRVE
while ePrevBreak<nEpsilonRVE
%% Pre-calculation of parameters for this step in the RVE run
    % Ensure crack locations are ordered and within x=[0,Lf]
    xDiscf(xDiscf>Lf)=xDiscf(xDiscf>Lf)-Lf;
    xDiscf=sort(xDiscf,3);

    % Updates interaction thickness of failed interactions (status=2)
    tv(statusv==2)=Inf;
    th(statush==2)=Inf;

    %Matrix of fibre thicknesses for interactions with k=1:4 neighbours
    Tf = aux_defineTfI(Nfi,Nfj,Df,tv,th);

    %Mapping fibre matrices to Interactions matrices
    [Tav,Tbv,Tah,Tbh] = aux_mapping_Tf_to_TI(Tf);
    
    statusv(Tav==0)=2;
    statusv(Tbv==0)=2;
    statush(Tah==0)=2;
    statush(Tbh==0)=2;
    
    %Mapping fibre discontinuities  3D matrices to Interactions 3D matrices
    [xDiscav,xDiscbv,xDiscah,xDiscbh] = aux_mapping_f_to_I(xDiscf);
    
    %Vector of strains which need updating (from ePrevBreak onwards)
    EpsilonRVENow=EpsilonRVE(ePrevBreak+1:end);
    nEpsilonRVENow=length(EpsilonRVENow);
        
    %% Calculates stress-strain curves for interactions which need update 
    %%(outputs interaction stresses interpolated to EpsilonRVEGrid)

    % Runs all vertical interactions (only returns for e>ePrevBreak)
    [SigmavRVE(:,:,ePrevBreak+1:end),statusv] = Interactions_SS(...
        xDiscav,xDiscbv,...
        SumEv,DeltaEv,tv,Gv,Sv,Giicv,...
        Lf,Eav,Ebv,Tav,Tbv,EpsilonRVENow,SigmavRVE(:,:,ePrevBreak+1:end),...
        statusv,EpsilonRVE(2));

    % Runs all horizontal interactions (only returns for e>ePrevBreak)
    [SigmahRVE(:,:,ePrevBreak+1:end),statush] = Interactions_SS(...
        xDiscah,xDiscbh,...
        SumEh,DeltaEh,th,Gh,Sh,Giich,...
        Lf,Eah,Ebh,Tah,Tbh,EpsilonRVENow,SigmahRVE(:,:,ePrevBreak+1:end),...
        statush,EpsilonRVE(2));
   
    % Re-sets the status of all interactions to be updated back to 0
    statusv(statusv==1)=0;
    statush(statush==1)=0;
    
    %% Calculates fibre stresses (average, failed fibres & crit neighbours)
    % Map SSI from interactions to fibres (for e>ePrevBreak)

    [SigmaFavg(:,:,ePrevBreak+1:nEpsilonRVE),...
        kfCrit,DeNextBreak,iFBreak,jFBreak] = ...
        Fibres_SS(Nfi,Nfj,Df,Tf,Ef,Eneighbours,...
        SigmavRVE(:,:,ePrevBreak+1:end),SigmahRVE(:,:,ePrevBreak+1:end),...
        XfNow,nEpsilonRVENow);
    
    %% Adds fibre discontinuity and updates status of interactions
    
    for q=1:length(iFBreak)
        [statusv,statush,xDiscf(iFBreak(q),jFBreak(q),1:5),...
            XfNow(iFBreak(q),jFBreak(q))] = Fibre_Breaks(...
            Nfi,Nfj,Lf,iFBreak(q),jFBreak(q),xDiscf,...
            kfCrit(iFBreak(q),jFBreak(q),DeNextBreak),...
            Xf(iFBreak(q),jFBreak(q),:),statusv,statush);
    end
        
    %% Prepares for next step in EpsilonRVE with updated fibre breaks
    
    ePrevBreak=DeNextBreak+ePrevBreak+1;
    
end

%Calculates the overall stress-strain curve of the RVE
SigmaRVE=permute(sum(sum(SigmaFavg.*Df.^2,1),2),[3 1 2])/sum(sum(Df.^2,1),2)*Vf;

%% Fracture criterion

[iCluster,jCluster,eMax,XRVEMax] = ...
    SFCfracture(nfInCluster,SigmaFavg,Df,SigmaRVE,EpsilonRVE,Vf,Nfi,Nfj,nEpsilonRVE);

end
