function [x,DX,S,g,...
    nMax,XeRemoteAll,PropertiesAll,XePointsAll,ijActive,DXL,gCrit] = ...
    SPSLLoop(N,type,Gm,y,gm,Sm,Lpz,Lall,T,tm,Eb,Nx,Ng)
% Loops calculations for all values of L provided

% Initiating matrices
nLall=length(Lall);
XeRemoteAll=NaN(1,2*nLall); % eRemote(1),XRemote(1),...
PropertiesAll=NaN(nLall,5); % L(1),E0(1),Xmax(1),eatXmax(1),UOverlap(1);...
XePointsAll=NaN(nLall,6); % eXnL(1),eXsoft(1),eXcrack(1);... 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Loop for all L required %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for l=1:nLall
    L=Lall(l);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main calculations for Shear-Lag Response %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x,DX,S,g,nMax,colL,ijActive,gCrit] = ...
        SPSLResponse(N,type,Gm,y,gm,Sm,Lpz,L,T,tm,Eb,Nx,Ng);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Post-processing results (all resulting strains are in percentage) %%%
    [XeRemote,Properties,XePoints,DXL] =...
        SPPostProcProperties(N,L,T,tm,Eb,Ng,DX,g,colL,ijActive);

    % Update nr of lines in XeRemote
    XeRemoteAll(size(XeRemoteAll,1)+1:size(XeRemote,1),:)=NaN;
    
    XeRemoteAll(:,2*l-1:2*l)=XeRemote;
    PropertiesAll(l,:)=[L,Properties];
    XePointsAll(l,:)=XePoints;
end

end