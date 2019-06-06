function [XeRemote,Properties,XePoints,DXL] =...
    SPPostProcProperties(N,L,T,tm,Eb,Ng,DX,g,colL,ijActive)
% Posto-processing the main properties (suitable for several L)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall stress-strain curve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DXL=DX(sub2ind(size(DX),1:size(DX,1),colL'))';   % Fields at the
gL=g(sub2ind(size(g),1:size(g,1),colL'))';       % edge of the overlap 

XRemote=DXL/2;                      % Eq.7a
eRemote=DXL/(2*Eb)+tm*gL/(2*L);     % Eq.7a

XeRemote=[eRemote*100,XRemote];     % Stress vs. strain (percentage) curve

XeRemote(end+1,:)=XeRemote(end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of main properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Strength and strain at maximum stress
[Xmax,lineXmax]=max(XRemote);
eatXmax=XeRemote(lineXmax,1);

% Elastic stiffness
E0=XeRemote(2,2)/(XeRemote(2,1)/100);

% Energy dissipated by the overlap (within the model length 2L)
% (crack opening is tm*gL; fracture toughness is trapz(tm*gL,XRemote) )
UOverlap=trapz(tm*gL,XRemote)*2*T;

% Summary of properties
Properties=[E0,Xmax,eatXmax,UOverlap];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of particular points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XePoints=[XeRemote(Ng+2,:),...                  % Onset of non-linearity
    XeRemote(find(ijActive(:,2)==N,1),:),...    % Onset of softening
    XeRemote(find(ijActive(:,2)==N+1,1),:)];    % Fully formed crack tip

end