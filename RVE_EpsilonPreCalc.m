function     [factorSigmaItoFavg,factorSigmaItoFpeak,...
    SumTv,DeltaTv,SumTh,DeltaTh,...
    xDiscav,xDiscbv,xDiscah,xDiscbh,XfNow] = ...
    RVEEpsilonPreCalc(Nfi,Nfj,Df,Ef,Xf,Eneighbours,tv,th,xDiscf,nDiscf)
%% Fibre thicknesses & parameters for interactions with neighbours

%Matrix of fibre thicknesses for interactions with k=1:4 neighbours
Tf = defineTfI(Nfi,Nfj,Df,tv,th);

%Mapping fibre matrices to Interaction matrices
[Tav,Tbv,Tah,Tbh] = mapping_Tf_to_TI(Tf);

%Calculating Sums and Deltas of Thicknesses
[SumTv,DeltaTv] = SumDelta(Tav,Tbv);
[SumTh,DeltaTh] = SumDelta(Tah,Tbh);

%% Generates Tf*Df
% Normalised area (Tf*Df) for each fibre, for the 4 interactions
TtimesDf=Tf.*Df;

% Normalised area (Tf*Df) for the neighbours of each fibre
TtimesDneighbours=zeros(Nfi,Nfj,4);
%1, above:
TtimesDneighbours(2:end,:,1)=Tav(:,:).*Df(1:end-1,:);
%2, on the right:
TtimesDneighbours(:,1:end-1,2)=Tbh(:,:).*Df(:,2:end);
%3, below:
TtimesDneighbours(1:end-1,:,3)=Tbv(:,:).*Df(2:end,:);
%4, on the left:
TtimesDneighbours(:,2:end,4)=Tah(:,:).*Df(:,1:end-1);

% Calculates factor to go from SigmaI to SigmaF (Eq 18 and 19 in paper)
factorSigmaItoFavg=Ef.*(TtimesDf.^2+TtimesDneighbours.^2)./...
    (Ef.*TtimesDf.^2+Eneighbours.*TtimesDneighbours.^2);
factorSigmaItoFpeak=(TtimesDf.^2+TtimesDneighbours.^2)./(TtimesDf.^2);

%% Updates array of discontinuities and fibre strengths
% Mapping location of fibre discontinuities
[xDiscav,xDiscbv,xDiscah,xDiscbh] = mapping_f_to_I(xDiscf);

% Fibre strength sets as the minimum value of the 4 "fragments-to-be"
XfNow=Xf(:,:,nDiscf);

end



