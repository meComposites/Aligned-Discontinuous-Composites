function [li,DXLi]=SPLengthSing(typei,yi,g0i,gLi,S0i,SLi,DX0i,T,tm,Eb)
% Calculates the length of individual subdomains to develop shear strains
% between g0i (at x=x0i) to gLi (at x=xLi).
% Based on Tables 2 and 3.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating subdomain length for each type of subdomain %%%%%%%%%%%%%%%%%
switch typei
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Positive stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        li = 2/yi*atanh(yi*T./(2*(S0i+SLi)).*(sqrt(...
            DX0i.^2+(2/(yi*T))^2*(SLi.^2-S0i.^2))-DX0i));   %Table 2, Eq.a
    
        DXLi = sqrt(DX0i.^2+(2/(yi*T))^2*(SLi.^2-S0i.^2));  %Table 3, Eq.a
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Zero stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0
        li = T./(2*S0i).*...
            (sqrt(DX0i.^2+4*tm*Eb/T.*S0i.*(gLi-g0i))-DX0i); %Table 2, Eq.d
        
        DXLi = sqrt(DX0i.^2+4*tm*Eb/T.*S0i.*(gLi-g0i));     %Table 3, Eq.b
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Negative stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case -1
        li = 2/yi*atan(yi*T./(2*(S0i+SLi)).*...             %Table 2, Eq.f
            (sqrt(DX0i.^2+(2/(yi*T))^2*(S0i.^2-SLi.^2))-DX0i));

        DXLi = sqrt(DX0i.^2+(2/(yi*T))^2*(S0i.^2-SLi.^2));  %Table 3, Eq.c
       
end

end