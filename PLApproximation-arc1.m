function [] = PLApproximation( L )
%JOEL Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Recuperation of results via SL Inputs & PreCalc %%%%%%%%%%%%%%%%%%%%%

figure

for Lall=L
    [XeRemoteAll,~,~,~,fieldx,~,fieldS,~,~,gCrit]=SL(Lall);
    [gmIn,SmIn,Giicm,T,tm,Eb,Nx,Ng] = Inputs();
    [N,type,Gm,y,gm,Sm,Sgm] = PreCalc(gmIn,SmIn,Giicm,T,tm,Eb);
    Lpz = LpzMatrix(N,type,y,gm,Sm,T,tm,Eb)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Use a piecewise linear approximation at critical points %%%%%%%%%%%%%

    XeRemoteAllPL=XeRemoteAll(1:Ng+1:size(XeRemoteAll,1),:);
    
    hold on
    plot(XeRemoteAll(:,1),XeRemoteAll(:,2),'r')
    hold on
    plot(XeRemoteAllPL(:,1),XeRemoteAllPL(:,2))
    for i=1:size(XeRemoteAllPL,1)
        hold on
        plot(XeRemoteAllPL(i,1),XeRemoteAllPL(i,2),'*')
    end
    
end
end

