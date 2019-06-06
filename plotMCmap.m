function plotMCmap(data,Nfi)
avg=nanmean(data(:));
std=nanstd(data(:));

nMC=size(data,3);

dataSymm=nan(size(data,1),size(data,2),nMC*16);
dataSymm(:,:,1:nMC)=data;
dataSymm(:,:,1*nMC+1:2*nMC)=rot90(data,1);
dataSymm(:,:,2*nMC+1:3*nMC)=rot90(data,2);
dataSymm(:,:,3*nMC+1:4*nMC)=rot90(data,3);
dataSymm(:,:,4*nMC+1:8*nMC)=flip(dataSymm(:,:,1:4*nMC),1);
dataSymm(:,:,8*nMC+1:12*nMC)=flip(dataSymm(:,:,1:4*nMC),2);
dataSymm(:,:,12*nMC+1:16*nMC)=flip(dataSymm(:,:,4*nMC+1:8*nMC),2);
dataAvg=nanmean(dataSymm,3);

% Create figure
figure('InnerPosition',[100 100 700 700]);

% Create axes
axes1 = axes;
axis off
hold(axes1,'on');

nIntervals=100;
blueRedMap=ones(2*nIntervals,3);
blueRedMap(1:nIntervals,2:3)=repmat(linspace(0,1,nIntervals)',1,2);
blueRedMap(nIntervals+1:end,1:2)=repmat(linspace(1,0,nIntervals)',1,2);
colormap(blueRedMap);    

% Create surface
surface(dataAvg,'AlignVertexCenters','on','LineStyle','none');

% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes1,[Nfi/2 3/2*Nfi]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes1,[Nfi/2 3/2*Nfi]);
%box(axes1,'on');
% Set the remaining axes properties

set(axes1,'CLim',[avg-std avg+std]);
set(axes1,'Positio',[0 0 1 1]);

