function [statusvNew,statushNew,xDiscfF,XfNowij] = Fibre_Breaks(...
    Nfi,Nfj,Lf,iF,jF,xDiscf,kfCritij,Xfij,statusv,statush)
% Performs ... for broken fibre
% and returns the failed fibres and critical neighbours

%% Defines critical neighbour and interactions to be updated 
% (all interacting with fibres a and b in broken interaction)
statusvNew=statusv;
statushNew=statush;

switch kfCritij
    case 1 
        %Neighbour above
        iN=iF-1;    jN=jF;
        %V: above N to below F 
        statusvNew(max(iN-1,1):min(iF,Nfi-1),jF)=1;
        %H: left and right of N & F:
        statushNew(iN:iF,max(jF-1,1):min(jF,Nfj-1))=1;
        %Broken interaction: above F
        statusvNew(iF-1,jF)=2;

    case 3 
        %Neighbour below
        iN=iF+1;    jN=jF;
        %V: above F to below N 
        statusvNew(max(iF-1,1):min(iN,Nfi-1),jF)=1;
        %H: left and right of F & N:
        statushNew(iF:iN,max(jF-1,1):min(jF,Nfj-1))=1;
        %Broken interaction: below F
        statusvNew(iF,jF)=2;

    case 4 
        %Neighbour left
        iN=iF;    jN=jF-1;
        %H: left N to right F 
        statushNew(iF,max(jN-1,1):min(jF,Nfj-1))=1;
        %H: above and below of N & F:
        statusvNew(max(iF-1,1):min(iF,Nfi-1),jN:jF)=1;
        %Broken interaction: left of F
        statushNew(iF,jF-1)=2;

    case 2 
        %Neighbour right
        iN=iF;    jN=jF+1;
        %H: left F to right N 
        statushNew(iF,max(jF-1,1):min(jN,Nfj-1))=1;
        %H: above and below of F & N:
        statusvNew(max(iF-1,1):min(iF,Nfi-1),jF:jN)=1;
        %Broken interaction: right of F
        statushNew(iF,jF)=2;
end

%Ensures that no previously-broken interaction (=2) will be updated 
statusvNew(statusv==2)=2;
statushNew(statush==2)=2;

%% Adds fibre break

%segmentsI: segments in the (critical) interaction between F and N
%xDiscI: vector of ordered discontinuities in both fibres
%1st row: length of the segments / location of discontinuities
%2nd row: type of segment (0=broken F (a); 1=broken N (b); 0.5=SL) 
[segmentsI,xDiscI] = Interaction_segments(Lf,...
    xDiscf(iF,jF,:),xDiscf(iN,jN,:));

% Zero-ing length of BRoken (i.e. not Shear-Lag) segments
segmentsI(segmentsI(:,2)~=0.5,1)=0;
[~,critSL]=max(segmentsI(:,1));

% Gets x of the break in the critical neighbouring fibre SL segment
% (Neighbouring fibre has ID=1, broken Fibre has ID=0)
xDiscfF=xDiscf(iF,jF,:);
xDiscfF(1,1,5)=sum(...
    xDiscI(critSL:critSL+1,1).*...
    xDiscI(critSL:critSL+1,2));
xDiscfF=sort(xDiscfF,3);

XfNowij=Xfij(1,1,sum(~isnan(xDiscfF),3));

end