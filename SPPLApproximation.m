function [] = SPPLApproximation( L )
%JOEL Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Recuperation of results via SL Inputs & PreCalc %%%%%%%%%%%%%%%%%%%%%

ErrX=[];
ErrY=[];
for Lall=L
    figure
    [XeRemoteAll,~,~,~,fieldx,~,fieldS,~,~,gCrit]=SPSL(Lall);
    [gmIn,SmIn,Giicm,T,tm,Eb,Nx,Ng] = SPInputs();
    [N,type,Gm,y,gm,Sm,Sgm] = SPPreCalc(gmIn,SmIn,Giicm,T,tm,Eb);
    Lpz = SPLpzMatrix(N,type,y,gm,Sm,T,tm,Eb)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Use a piecewise linear approximation at critical points %%%%%%%%%%%%%
    
    XeRemoteAllPL=XeRemoteAll(1:Ng+1:size(XeRemoteAll,1),:);
    Nbr=(size(XeRemoteAll)+1)/Ng-1;
    for i=1:Nbr
        XeRemoteAll=XeRemoteAll([1:i*Ng,i*Ng+2:end],:);
    end
    hold on
    plot(XeRemoteAll(:,1),XeRemoteAll(:,2),'r')
    hold on
    plot(XeRemoteAllPL(:,1),XeRemoteAllPL(:,2))
    for i=1:size(XeRemoteAllPL,1)
        hold on
        plot(XeRemoteAllPL(i,1),XeRemoteAllPL(i,2),'*')
    end
    
    x=1:1:length(XeRemoteAll(:,1));
    XeRemoteAllPLErr(:,2)=XeRemoteAllPL(:,2);
    XeRemoteAllPLErr(:,1)=x(1:Ng:size(x,2));
    
    interp=interp1(XeRemoteAllPLErr(:,1),XeRemoteAllPL(:,2),x)';
    figure
    plot(x,abs((interp-XeRemoteAll(:,2))./XeRemoteAll(:,2)*100))
    
    figure
    plot(x,interp,x,XeRemoteAll(:,2),'*')
    
    ErrX=[ErrX Lall];
    ErrY=[ErrY max(abs((interp-XeRemoteAll(:,2))./XeRemoteAll(:,2)*100))];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Use a polynomial approximation at critical points %%%%%%%%%%%%%%%%%%%
    
    %     N=2
    %     Pol=polyfit(XeRemoteAllPL(:,1),XeRemoteAllPL(:,2),N)
    %     x=0:step:eInf;
    %     y=0;
    %     for i=0:N
    %         y=y+Pol(i+1)*x.^(N-i);
    %     end
    %     hold on
    %     plot(x,y)
    
end

if length(L)>1
    figure
    plot(ErrX,ErrY)
end
end

