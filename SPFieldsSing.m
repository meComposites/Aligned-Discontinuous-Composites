function [xi,DX,S,g]=SPFieldsSing(typei,Gmi,yi,g0i,S0i,DX0i,li,T,tm,Eb,Nx,Ng)
% Calculates stress and strain fields 
% for a single subdomain (Table 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating location of x-points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xi=SPlinspace2(0,li,Nx+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating matrices with fields at xi=0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g0i=repmat(g0i,1,Nx+1);
DX0i=repmat(DX0i,1,Nx+1);
S0i=repmat(S0i,1,Nx+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating fields for each type of subdomain %%%%%%%%%%%%%%%%%%%%%%%%%%%
switch typei
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Positive stiffness, Table 1 Eq.b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1      
        DX=DX0i.*cosh(yi*xi)+2/(yi*T)*S0i.*sinh(yi*xi);
        S=S0i.*cosh(yi*xi)+yi*T/2*DX0i.*sinh(yi*xi);
        g=g0i+S0i/Gmi.*(cosh(yi*xi)-1)+DX0i/(yi*tm*Eb).*sinh(yi*xi);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Zero stiffness, Table 1 Eq.d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0      
        DX=DX0i+2/T*S0i.*xi;
        S=S0i.*ones(Ng+1,Nx+1);
        g=g0i+S0i./(T*tm*Eb).*xi.^2+DX0i/(tm*Eb).*xi;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Negative stiffness, Table 1 Eq.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case -1     
        DX=DX0i.*cos(yi*xi)+2/(yi*T)*S0i.*sin(yi*xi);
        S=S0i.*cos(yi*xi)-yi*T/2*DX0i.*sin(yi*xi);
        g=g0i+S0i/Gmi.*(1-cos(yi*xi))+DX0i/(yi*tm*Eb).*sin(yi*xi);
end

end