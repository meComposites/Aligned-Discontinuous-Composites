function [g0Crit] = SPTransSing(typei,Gmi,yi,gmi,Smi,L,T,tm,Eb)
% Calculates the transition strain at the centre of the overlap 
% (only one subdomain in the overlapping length, hence transition 
% corresponds to the activation of the next subdomain).
% Based on Table 4, Eqs.b-d.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating transition strains for each type of subdomain %%%%%%%%%%%%%%%
switch typei
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Positive stiffness, Table 4 Eq.b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        g0Crit=gmi-Smi/Gmi+Smi/(Gmi*cosh(yi*L));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Zero stiffness, Table 4 Eq.c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 0
        g0Crit=gmi-Smi*L^2/(T*tm*Eb);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Negative stiffness, Table 4 Eq.d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case -1
        g0Crit=gmi-Smi/(Gmi*cos(yi*L))+Smi/Gmi;
        
end

end