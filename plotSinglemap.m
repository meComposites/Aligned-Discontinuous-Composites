function plotSinglemap(data,Nfi,iCluster,jCluster,nfInCluster)
avg=nanmean(data(:));
std=nanstd(data(:));

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
surface(data,'AlignVertexCenters','on','LineStyle','none',...
    'MarkerFaceColor','flat',...
    'MarkerSize',10,...
    'Marker','o',...
    'LineStyle','none',...
    'FaceColor','none')


% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes1,[0.5 Nfi+0.5]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes1,[0.5 Nfi+0.5]);
%box(axes1,'on');
% Set the remaining axes properties

set(axes1,'CLim',[avg-3*std avg+3*std]);
set(axes1,'Positio',[0 0 1 1]);

hold on

xCluster=[jCluster-0.5;...
    jCluster-0.5;...
    jCluster+nfInCluster-0.5;...
    jCluster+nfInCluster-0.5;...
    jCluster-0.5];
yCluster=[iCluster-0.5;...
    iCluster+nfInCluster-0.5;...
    iCluster+nfInCluster-0.5;...
    iCluster-0.5;...
    iCluster-0.5];

%plot(xCluster,yCluster,'LineWidth',5.0,'Color',[0 0 0]);

hold off

