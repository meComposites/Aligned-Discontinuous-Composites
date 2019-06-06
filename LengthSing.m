function [li,Pi]=LengthSing(typei,yi,g0i,gLi,S0i,SLi,P0i,ET,DT,tm,EE,DE)
% Calculates the length of individual subdomains to develop shear strains
% between g0i (at x=x0i) to gLi (at x=xLi).
% Based on Tables 2 and 3.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating subdomain length for each type of subdomain %%%%%%%%%%%%%%%%%
switch typei
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Positive stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        T=((ET^2-DT^2)/(2*ET));
        li = 2/yi*atanh(yi*T./(2*(S0i+SLi)).*(sqrt(...
            P0i.^2+(2/(yi*T))^2*(SLi.^2-S0i.^2))-P0i));   %Table 2, Eq.a
        
        li(S0i<2*eps(SLi))=1/yi*log(2*SLi./S0i(S0i<2*eps(SLi)));
        
        Pi = sqrt(P0i.^2+(2/(yi*T))^2*(SLi.^2-S0i.^2));  %Table 3, Eq.a
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Zero stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0
        T=((ET^2-DT^2)/(2*ET));
        Eb=2*ET*(EE-DE)/2*(EE+DE)/2/(EE*ET+DE*DT);
        %Eb=ET*(EE^2-DE^2)/((ET-DT)*(EE-DE)+(ET+DT)*(EE+DE));
        li = T./(2*S0i).*...
            (sqrt(P0i.^2+4*tm*Eb/T.*S0i.*(gLi-g0i))-P0i); %Table 2, Eq.d
        
        Pi = sqrt(P0i.^2+4*tm*Eb/T.*S0i.*(gLi-g0i));     %Table 3, Eq.b
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Negative stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case -1
        T=((ET^2-DT^2)/(2*ET));
        li = 2/yi*atan(yi*T./(2*(S0i+SLi)).*(sqrt(P0i.^2+(2/(yi*T))^2*(S0i.^2-SLi.^2))-P0i)); %Table 2, Eq.f
        li(S0i==0)=pi/(2*yi);
        
%         if S0i==0
%            li=pi/(2*yi); 
%            
%         else
%            li = 2/yi*atan(yi*T./(2*(S0i+SLi)).*(sqrt(P0i.^2+(2/(yi*T))^2*(S0i.^2-SLi.^2))-P0i)); %Table 2, Eq.f
%         end
        
        Pi = sqrt(P0i.^2+(2/(yi*T))^2*(S0i.^2-SLi.^2));  %Table 3, Eq.c
       
end

end