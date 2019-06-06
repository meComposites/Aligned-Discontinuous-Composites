function [] = linkEpsSig( L )
%JOEL Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Recuperation of results via SL Inputs & PreCalc %%%%%%%%%%%%%%%%%%%%%
figure
for Lall=L
    [XeRemoteAll,~,~,~,~,~,~,~,~,~]=SL(Lall);
    [gmIn,SmIn,Giicm,T,tm,Eb,~,Ng] = Inputs();
    [N,type,Gm,y,gm,Sm,Sgm] = PreCalc(gmIn,SmIn,Giicm,T,tm,Eb);
    Lpz = LpzMatrix(N,type,y,gm,Sm,T,tm,Eb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Find theoritical function between epsilon and sigma sk=[1] %%%%%%%%%%
        G=Gm(1);
        alpha=1/(y(1)*Lall*tanh(y(1)*Lall));
        Emat=Eb/(1+alpha);

        EpsEndLin=XeRemoteAll(Ng,1); % end of linear part
        SigEndLin=EpsEndLin*Emat/100;
        eTh=[0 EpsEndLin];
        sTh=[0 SigEndLin];
    
        hold on
        plot(XeRemoteAll(:,1),XeRemoteAll(:,2),'r')
        hold on
        plot(EpsEndLin,SigEndLin,'*')
        hold on
        plot(eTh,sTh)
        axis([0,3*EpsEndLin,0,3*SigEndLin])
        
       
end
end

