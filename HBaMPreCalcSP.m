function [SumT,DeltaT,SumE,DeltaE,gI1,gI2,Y1,Y2,I1_lt,Teq,eta,K] = HBaMPreCalc(...
    Ta,Ea,Tb,Eb,SI,GiicI,tI,GI,opts)

signDeltaK=(Tb.*Eb-Ta.*Ea)./abs((Tb.*Eb-Ta.*Ea));

% Sum and difference of Tf and Ef
SumT=Ta+Tb;
DeltaT=signDeltaK.*(Tb-Ta);

SumE=Ea+Eb;
DeltaE=signDeltaK.*(Eb-Ea);

K=2*(DeltaE.*SumT)./(SumE.*SumT+DeltaE.*DeltaT); % kappa
eta=(2*SumT./(SumT-DeltaT)+K)./(2*SumT./(SumT+DeltaT)-K); % eta=(EB.TE)/(EA.TA)

Teq=((SumT.^2-DeltaT.^2)./(2*SumT)); % Equivalent thickness

%Matrix constitutive law
gI1=SI./GI;
gI2=2*GiicI./(tI.*SI);
GI2=SI./(gI1-gI2);

Y1=sqrt(8*GI./tI.*(SumE.*SumT+DeltaE.*DeltaT)./...
    ((SumE.^2-DeltaE.^2).*(SumT.^2-DeltaT.^2)));
Y2=sqrt(-8*GI2./tI.*(SumE.*SumT+DeltaE.*DeltaT)./...
    ((SumE.^2-DeltaE.^2).*(SumT.^2-DeltaT.^2)));

% Other calculations SP - I DONT KNOW WHAT THEY ARE!!
I1_PsiL=2*SI./(Y2.*Teq);
I1_PsiR=-I1_PsiL./eta;
I1_tR=sqrt(SI.^2-I1_PsiR.^2.*Y2.^2.*Teq.^2/4);
I1_lt=2./Y2.*(atan(1)+atan(sqrt((SI-I1_tR)./(SI+I1_tR))));


%[P2_g0_212,~,~]= fminbnd(@P2_F_212,0,gI1,opts);

end

