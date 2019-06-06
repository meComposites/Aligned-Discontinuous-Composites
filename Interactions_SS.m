function [SigmaIAtRVEGridNow,statusI] = Interactions_SS(xDiscaI,xDiscbI,...
    EE,DE,tI,GI,SI,GiicI,...
    Lf,Ea,Eb,Ta,Tb,EpsilonRVEGridNow,SigmaIAtRVEGridNow,statusI,DeltaEpsilonRVE)
%Calculates the stress-strain curve of the interactions which need updating

%% Calculates fibre thicknesses & parameters for interactions

%Calculating Sums and Deltas of Thicknesses
[ET,DT] = aux_SumDelta(Ta,Tb);

nDiscaI=sum(~isnan(xDiscaI),3);
nDiscbI=sum(~isnan(xDiscbI),3);

%% For all interactions which need updating:

[iItoUpdate,jItoUpdate] = find(statusI==1);

for q=1:length(iItoUpdate)
    i=iItoUpdate(q);
    j=jItoUpdate(q);
    
    %% Creates length of segments in the interaction
    [segmentsI,~] = Interaction_segments(Lf,...
        xDiscaI(i,j,1:nDiscaI(i,j)),...
        xDiscbI(i,j,1:nDiscbI(i,j)));
    
    % Vector with SL segments (average=0.5)
    lSLI=sort(segmentsI(segmentsI(:,2)==0.5,1));

    % Total length of Br segment with both ends in the "a" fibre
    lBra=sum(segmentsI(segmentsI(:,2)==0,1));
    % Total length of Br segment with both ends in the "b" fibre
    lBrb=sum(segmentsI(segmentsI(:,2)==1,1));

    %% Runs  SS_HBaM for all overlaps Lo in this interaction,
    %  if the interaction is not broken
    if lSLI(1)>0.01
    
        nSLI=length(lSLI);
        ESLI=zeros(1,nSLI);
        deltaeSLI=ESLI;
        EpsilonSLI=nan(50,nSLI);
        SigmaSLI=EpsilonSLI;
        for k=1:nSLI
            signDeltaK=abs(Eb(i,j)*Tb(i,j)-Ea(i,j)*Ta(i,j))/...
                (Eb(i,j)*Tb(i,j)-Ea(i,j)*Ta(i,j));
            
%              % Calculates overlap - old method
%              [ESLI(k),deltaeSLI(k),EpsilonNow,SigmaNow] = ...
%                  SS_HBaM_BiLin(...
%                  EE(i,j),signDeltaK*DE(i,j),ET(i,j),signDeltaK*DT(i,j),tI(i,j),...
%                  lSLI(k),SI(i,j),GiicI(i,j),GI(i,j));
            
            %Calculates overlap - new method
            [ESLI(k),deltaeSLI(k),EpsilonNow,SigmaNow] = ...
               HBaMBiLinSP(...
               EE(i,j),signDeltaK*DE(i,j),ET(i,j),signDeltaK*DT(i,j),tI(i,j),...
               lSLI(k),SI(i,j),GiicI(i,j),GI(i,j),500);    
            
            %plot(EpsilonNow,SigmaNow,EpsilonNow2,SigmaNow2);
            
            EpsilonSLI(1:length(SigmaNow),k)=EpsilonNow;
            SigmaSLI(1:length(SigmaNow),k)=SigmaNow;
        end        
        %% Combines all segments of the interaction in series, 
        %  and outputs to EpsilonGrid
        SigmaIAtRVEGridNow(i,j,:) = combineSeriesI(...
            nSLI,lSLI,EpsilonSLI,SigmaSLI,ESLI,deltaeSLI,...
            lBra,lBrb,Ea(i,j),Eb(i,j),Lf,ET(i,j),Ta(i,j),Tb(i,j),...
            EpsilonRVEGridNow,DeltaEpsilonRVE);
    else
        statusI(i,j)=2;
    end

end

%% For all interactions which need have failed:

[iIFailed,jIFailed] = find(statusI==2);

for q=1:length(iIFailed)
    i=iIFailed(q);
    j=jIFailed(q);
    SigmaIAtRVEGridNow(i,j,:)=0;
end

end
