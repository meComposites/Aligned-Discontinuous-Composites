function [gi,Si,Pi]=field_2Stiffness_nonIterative(State,k_m,EE,DE,G,Sm,gm,Y,K,ET,DT,tm,g0,S0,P0,Xinf,Dli,typei)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating fields for each type of subdomain %%%%%%%%%%%%%%%%%%%%%%%%%%%
switch typei
    
    case 1
%         DXi=DX0*cosh(Y(i_f,j_f,k_m)*(Dli))+2/(Y(i_f,j_f,k_m)*T)*S0*sinh(Y(i_f,j_f,k_m)*(Dli))+(DXL*DEij(j_f,i_f)-X0i(j_f)*E(i_f)+X0i(i_f)*E(j_f))/Eij(i_f,j_f)*(cosh(Y(i_f,j_f,k_m)*(Dli))-1);
%         gi=g0+Eij(i_f,j_f)/(tm*E(i_f)*E(j_f))*(((DXL*DEij(j_f,i_f)-X0i(j_f)*E(i_f)+X0i(i_f)*E(j_f))/(Eij(i_f,j_f)*Y(i_f,j_f,k_m))+DX0/Y(i_f,j_f,k_m)).*sinh(Y(i_f,j_f,k_m)*(Dli))+(2*S0/((Y(i_f,j_f,k_m)^2)*T).*(cosh(Y(i_f,j_f,k_m)*(Dli))-1)));
%         Si=S0.*cosh(Y(i_f,j_f,k_m)*(Dli))+Y(i_f,j_f,k_m)*T/2*(DX0+(DXL*DEij(j_f,i_f)-X0i(j_f)*E(i_f)+X0i(i_f)*E(j_f))/(Eij(i_f,j_f))).*sinh(Y(i_f,j_f,k_m)*(Dli));
        
        Pi=P0*cosh(Y(k_m)*(Dli))+4*ET/(ET^2-DT^2)*S0/Y(k_m)*sinh(Y(k_m)*(Dli));
        Si=(ET^2-DT^2)/(4*ET)*P0*Y(k_m)*sinh(Y(k_m)*(Dli))+S0*cosh(Y(k_m)*(Dli));
        gi=g0+2/(tm*ET)*(EE*ET+DE*DT)/(EE^2-DE^2)*((P0)/(Y(k_m))*sinh(Y(k_m)*(Dli))+4*ET/(ET^2-DT^2)*S0/Y(k_m)^2*(cosh(Y(k_m)*(Dli))-1));
        
%         Pi2=P0*cosh(Y(k_m)*(-Dli))+4*ET/(ET^2-DT^2)*S0/Y(k_m)*sinh(Y(k_m)*(-Dli));
        
    case 0
%         DXi=DX0;
%         Si=Sm(end);
%         gi=g0+Dli/(E(i_f)*E(j_f)*tm)*(DX0*Eij(i_f,j_f)+DXL*DEij(j_f,i_f)+E(j_f)*X0i(i_f)-E(i_f)*X0i(j_f));

        if S0 > 0   % Zero stiffness
            Pi=P0+4*ET/(ET^2-DT^2)*S0*Dli;
            Si=S0;
            gi=g0+2/(tm*ET)*(EE*ET+DE*DT)/(EE^2-DE^2)*((P0)*Dli+2*ET/(ET^2-DT^2)*S0*Dli^2);
        else        % Fully debonded
            Pi=P0;
            Si=0;
            gi=g0+2/(tm*ET)*(EE*ET+DE*DT)/(EE^2-DE^2)*((P0)*Dli);
        end
        
    case -1
%         DXi=DX0*cos(Y(i_f,j_f,k_m)*(Dli))+2/(Y(i_f,j_f,k_m)*T)*S0*sin(Y(i_f,j_f,k_m)*(Dli))+(DXL*DEij(j_f,i_f)-X0i(j_f)*E(i_f)+X0i(i_f)*E(j_f))/Eij(i_f,j_f)*(cos(Y(i_f,j_f,k_m)*(Dli))-1);
%         gi=g0+Eij(i_f,j_f)/(tm*E(i_f)*E(j_f))*(((DXL*DEij(j_f,i_f)-X0i(j_f)*E(i_f)+X0i(i_f)*E(j_f))/(Eij(i_f,j_f)*Y(i_f,j_f,k_m))+DX0/Y(i_f,j_f,k_m)).*sin(Y(i_f,j_f,k_m)*(Dli))+(2*S0/((Y(i_f,j_f,k_m)^2)*T).*(1-cos(Y(i_f,j_f,k_m)*(Dli)))));
%         Si=S0.*cos(Y(i_f,j_f,k_m)*(Dli))-Y(i_f,j_f,k_m)*T/2*(DX0+(DXL*DEij(j_f,i_f)-X0i(j_f)*E(i_f)+X0i(i_f)*E(j_f))/(Eij(i_f,j_f))).*sin(Y(i_f,j_f,k_m)*(Dli));

        Pi=P0*cos(Y(k_m)*(Dli))+4*ET/(ET^2-DT^2)*S0/Y(k_m)*sin(Y(k_m)*(Dli));
        Si=-(ET^2-DT^2)/(4*ET)*P0*Y(k_m)*sin(Y(k_m)*(Dli))+S0*cos(Y(k_m)*(Dli));
        gi=g0+2/(tm*ET)*(EE*ET+DE*DT)/(EE^2-DE^2)*((P0)/(Y(k_m))*sin(Y(k_m)*(Dli))-4*ET/(ET^2-DT^2)*S0/Y(k_m)^2*(cos(Y(k_m)*(Dli))-1));
end


end