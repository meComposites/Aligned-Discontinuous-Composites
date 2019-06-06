function [fieldx,fieldg,fieldS,fieldXB] =...
    SPPostProcFields(Nx,nMax,x,DX,S,g,DXL)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full fields in the overlap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stresses in brick B
DXLrep=repmat(DXL,1,(nMax+1)*(Nx+1));
XBField=(DXLrep+DX)/2;      % Eq. 9d

% Symmetry conditions, Eq. 6b
fieldx=[-x(1:end,end:-1:1),x];
fieldS=[S(1:end,end:-1:1),S];
fieldg=100*[g(1:end,end:-1:1),g];
fieldXB=[DXLrep-XBField(1:end,end:-1:1),XBField]; 

end