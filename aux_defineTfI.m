%% Defines matrix of fibre thicknesses for the interaction with neighbours, TfI
function TfI = aux_defineTfI(Nfi,Nfj,Df,tv,th)

%Calculating distance from fibre(i,j) to its k=1:4 neighbours (3D matrices):
distfNeighbours=inf(Nfi,Nfj,4);
%1, above:
distfNeighbours(2:end,:,1)=Df(1:end-1,:)/2+tv;
%2, on the right:
distfNeighbours(:,1:end-1,2)=Df(:,2:end)/2+th;
%3, below:
distfNeighbours(1:end-1,:,3)=Df(2:end,:)/2+tv;
%4, on the left
distfNeighbours(:,2:end,4)=Df(:,1:end-1)/2+th;

%Define lambda (load sharing coefficient) for each fibre
y=1./sum(1./distfNeighbours.^2,3);

%Calculate fibre thicknesses for the interaction with each neighbour
TfI=repmat(y.*Df,1,1,4).*1./distfNeighbours.^2;

%Removes NaN from TfI due to fully broken fibre (from all neighbours)
TfI(isnan(TfI))=0;

end

