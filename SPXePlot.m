function [] = SPXePlot( Lall )
%JOEL Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Recuperation of results via SL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Results,~,XePointsAll]=SPSL(Lall);

% figure
% surface(fieldS,'LineStyle','none')
% colorbar
N=size(Results,2)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Definition and plot of results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ResultsAll=[];
AbsALL=[];

for i=0:1:N-1
    ResultsAll=[ResultsAll, Results(:,2*i+2)];
    AbsALL=[AbsALL, Results(:,2*i+1)];
end

figure;
plot(AbsALL,ResultsAll)
%plot(AbsALL,ResultsAll,'k')
title('Strain-stress curve')
xlabel('Strain (%)')
ylabel('Stress (Mpa)')
legend('L=0.5','L=1','L=10');


%%%%% Threshold points %%%%%
%for i=1:1:N
%    hold on
%    plot(XePointsAll(i,1),XePointsAll(i,2),'*',XePointsAll(i,3),XePointsAll(i,4),'*',XePointsAll(i,5),XePointsAll(i,6),'*')
%end

end

