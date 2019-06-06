function [li,DXLi]=LengthSingCentre(typei,yi,g0i,gLi,S0i,SLi,ET,DT,tm,EE,DE)
% Calculates the length of the central subdomain to develop shear strains
% between g0i (at x=x0i) to gLi (at x=xLi).
% Particular case of LengthSing, required due to numerical errors when
% L>Lpz and typei=1.
% Based on Tables 2 and 3.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating subdomain length for each type of central subdomain %%%%%%%%%
switch typei
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Positive stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        T=((ET^2-DT^2)/(2*ET));
        li=2/yi*atanh(sqrt((SLi-S0i)./(SLi+S0i)));          %Table 2, Eq.b
        
        li(S0i<2*eps(SLi))=1/yi*log(2*SLi./S0i(S0i<2*eps(SLi)));
        
        DXLi = 2/(yi*T)*sqrt(SLi.^2-S0i.^2);                %Table 3, Eq.a
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Zero stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0
        T=((ET^2-DT^2)/(2*ET));
        Eb=ET*(EE^2-DE^2)/((ET-DT)*(EE-DE)+(ET+DT)*(EE+DE));
        li = sqrt(T*tm*Eb.*(gLi-g0i)./S0i);                 %Table 2, Eq.d
        
        DXLi = sqrt(4*tm*Eb/T.*S0i.*(gLi-g0i));             %Table 3, Eq.b
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Negative stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case -1
        T=((ET^2-DT^2)/(2*ET));
        li = 2/yi*atan(sqrt((S0i-SLi)./(S0i+SLi)));         %Table 2, Eq.f
        
        li(S0i==0)=pi/(2*yi);                               %Table 2, Eq.h
        
        DXLi = 2/(yi*T)*sqrt(S0i.^2-SLi.^2);                %Table 3, Eq.c
       
end

end