% clear
% close all
dbstop if error
%% Inputs
%load('RVE_inputs.mat');
load('RVE_inputsICCM.mat');

%% Controlling the random generator
% rng(1);

Nfi=18; Nfj=387;

nfInCluster=2;
nMC=100;

Lf=3;

ncTarget=nc;

distNearestMC=zeros(Nfi,Nfj,nMC);
fibreStrengthMC=zeros(Nfi,Nfj,nMC);
fibreStiffnessMC=zeros(Nfi,Nfj,nMC);
SIminMC=zeros(Nfi,Nfj,nMC);
LominMC=zeros(Nfi,Nfj,nMC);

SigmaRVE=zeros(nEpsilonRVE,nMC);
eMax=zeros(1,nMC);
XMax=eMax; iCluster=eMax; jCluster=eMax;

for MC=1:nMC
    rng(MC);
    
    %% Calculating input parameters for RVE %%%%%
    % Random inputs for RVE microstructure
    [nc,xEndf,Typef,Df,Ef,Xf,...
        tv0,th0,Gv,Gh,Sv,Sh,Giicv,Giich] ...
        = RVE_generate...
        (Nfi,Nfj,Vf,Lf,ncTarget,...
        Dg,Dc,EgAvg,EgCoV,EcAvg,EcCoV,...
        XgAvg,XgCoV,XcAvg,XcCoV,...
        SIAvg,SICoV,GIAvg,GiicIAvg);
    
    %Ef=EcAvg.*(1+randn(size(Typef)).*0.001);
    %tv0=mean(tv0(:))*ones(Nfi-1,Nfj);
    %th0=mean(th0(:))*ones(Nfi,Nfj-1);
    %Sv=mean(Sv(:))*ones(Nfi-1,Nfj);
    %Sh=mean(Sh(:))*ones(Nfi,Nfj-1);
    %Xf(:,:,1:4)=sum(sum(Xf(:,:,1),1),2)/(Nfi*Nfj);
    %Xf(:,:,2)=sum(sum(Xf(:,:,2),1),2)/(Nfi*Nfj);
    %Xf(:,:,3)=sum(sum(Xf(:,:,3),1),2)/(Nfi*Nfj);
    %Xf(:,:,4)=sum(sum(Xf(:,:,4),1),2)/(Nfi*Nfj);
    
    [SigmaRVE(:,MC),eMax(1,MC),XMax(1,MC),iCluster(1,MC),jCluster(1,MC)]=...
        RVE_SS(nc,xEndf,Typef,Df,Ef,Xf,...
        tv0,th0,Gv,Gh,Sv,Sh,Giicv,Giich,...
        Lf,nEpsilonRVE,EpsilonRVEMax,opts,Vf,nfInCluster,Nfi,Nfj);
    
    %% Shear strength with neighbour
        SIfNeighbours=inf(Nfi,Nfj,4);
        %1, above:
        SIfNeighbours(2:end,:,1)=Sv;
        %2, on the right:
        SIfNeighbours(:,1:end-1,2)=Sh;
        %3, below:
        SIfNeighbours(1:end-1,:,3)=Sv;
        %4, on the left
        SIfNeighbours(:,2:end,4)=Sh;
        SImin=min(SIfNeighbours,[],3);   
        
           %% Overlap length with neighbour
        LoNeighbours=Lf/2*ones(Nfi,Nfj,4);
        %1, above:
        LoNeighbours(2:end,:,1)=xEndf(2:end,:)-xEndf(1:end-1,:);
        %2, on the right:
        LoNeighbours(:,1:end-1,2)=xEndf(:,1:end-1)-xEndf(:,2:end);
        %3, below:
        LoNeighbours(1:end-1,:,3)=xEndf(1:end-1,:)-xEndf(2:end,:);
        %4, on the left
        LoNeighbours(:,2:end,4)=xEndf(:,2:end)-xEndf(:,1:end-1);
        
        LoNeighbours=abs(LoNeighbours);
        LoNeighbours2=Lf-LoNeighbours;
        LoNeighbours=min(cat(4,LoNeighbours,LoNeighbours2),[],4);
        Lomin=min(LoNeighbours,[],3);
    
    %% Distance to nearest neighbour
        spacingfNeighbours=inf(Nfi,Nfj,4);
        %1, above:
        spacingfNeighbours(2:end,:,1)=tv0;
        %2, on the right:
        spacingfNeighbours(:,1:end-1,2)=th0;
        %3, below:
        spacingfNeighbours(1:end-1,:,3)=tv0;
        %4, on the left
        spacingfNeighbours(:,2:end,4)=th0;
        distNearest=min(spacingfNeighbours,[],3);

        distNearestMC(:,:,MC)=distNearest;
        
        
        fibreStrengthMC(:,:,MC)=Xf(:,:,1);
        fibreStiffnessMC(:,:,MC)=Ef;
        SIminMC(:,:,MC)=SImin;
        LominMC(:,:,MC)=Lomin;
        
MC

end

distNearestMaster=nan(Nfi*2-nfInCluster,2*Nfj-nfInCluster,nMC);
fibreStrengthMaster=nan(Nfi*2-nfInCluster,2*Nfj-nfInCluster,nMC);

SIminMaster=nan(Nfi*2-nfInCluster,2*Nfj-nfInCluster,nMC);

LominMaster=nan(Nfi*2-nfInCluster,2*Nfj-nfInCluster,nMC);


fibreStiffnessMaster=nan(Nfi*2-nfInCluster,2*Nfj-nfInCluster,nMC);

for q=1:nMC
   distNearestMaster((Nfi-nfInCluster+1:Nfi*2-nfInCluster)-iCluster(q)+1,...
       (Nfj-nfInCluster+1:Nfj*2-nfInCluster)-jCluster(q)+1,q)=distNearestMC(:,:,q); 
   
   fibreStiffnessMaster((Nfi-nfInCluster+1:Nfi*2-nfInCluster)-iCluster(q)+1,...
       (Nfj-nfInCluster+1:Nfj*2-nfInCluster)-jCluster(q)+1,q)=fibreStiffnessMC(:,:,q); 
   
   fibreStrengthMaster((Nfi-nfInCluster+1:Nfi*2-nfInCluster)-iCluster(q)+1,...
       (Nfj-nfInCluster+1:Nfj*2-nfInCluster)-jCluster(q)+1,q)=fibreStrengthMC(:,:,q); 
   
   SIminMaster((Nfi-nfInCluster+1:Nfi*2-nfInCluster)-iCluster(q)+1,...
       (Nfj-nfInCluster+1:Nfj*2-nfInCluster)-jCluster(q)+1,q)=SIminMC(:,:,q);
   
   LominMaster((Nfi-nfInCluster+1:Nfi*2-nfInCluster)-iCluster(q)+1,...
       (Nfj-nfInCluster+1:Nfj*2-nfInCluster)-jCluster(q)+1,q)=LominMC(:,:,q);
      
end

1;
