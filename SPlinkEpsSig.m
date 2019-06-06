function [] = SPlinkEpsSig( Lall )
%JOEL Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Recuperation of results via SL Inputs & PreCalc %%%%%%%%%%%%%%%%%%%%%
figure
for L=Lall
    [XeRemoteAll,PropertiesAll,~,~,~,~,~,~,~,~]=SPSL(L);
    [gmIn,SmIn,Giicm,T,tm,Eb,~,Ng] = SPInputs();
    [N,type,Gm,y,gm,Sm,Sgm] = SPPreCalc(gmIn,SmIn,Giicm,T,tm,Eb);
    Lpz = SPLpzMatrix(N,type,y,gm,Sm,T,tm,Eb);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Find theoritical function between epsilon and sigma sk=[1] %%%%%%%%%%
    
    alpha=1/(y(1)*L*tanh(y(1)*L));
    Emat=Eb/(1+alpha);
    
    EpsEndLin=XeRemoteAll(Ng+1,1); % end of linear part
    SigEndLin=EpsEndLin*Emat/100;
    eTh=[0 EpsEndLin];
    sTh=[0 SigEndLin];
    
    hold on
    plot(XeRemoteAll(:,1),XeRemoteAll(:,2),'-.k')
    hold on
    plot(EpsEndLin,SigEndLin,'*')
    hold on
    plot(eTh,sTh)
    %axis([0,3*EpsEndLin,0,3*SigEndLin])
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Find theoritical function between epsilon and sigma sk=[1,2] %%%%%%%%
    
    % Parameter: x1, varying from 0 or L-lpz[1] to L
    
    x1Min=max(0,L-Lpz(1,1));
    x1Max=L;
    x1StepNumber=100;
    x1Step=(x1Max-x1Min)/x1StepNumber;
    x1=x1Max:-x1Step:x1Min;
    tau1=Sgm(2,2);
    tau0=tau1./cosh(y(1).*x1);
    
    sigma12=tau0./(y(1)*T).*sinh(y(1).*x1)+tau1./T.*(L-x1);
    epsilon12=100*(1/(T*Eb).*(tau0./y(1).*sinh(y(1).*x1)+tau1.*(L-x1))+(tm.*tau0)./(2*L*Gm(1)).*(cosh(y(1).*x1))+tau1/(2*L*T*Eb).*(L-x1).^2+(L-x1).*tau0/(Eb*L*y(1)*T).*sinh(y(1).*x1));
    
    hold on
    plot(epsilon12,sigma12,'r')
    
    hold on
    plot(epsilon12(end),sigma12(end),'r*')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Find theoritical function between epsilon and sigma sk=[2] %%%%%%%%%%
    
    % Parameter: gamma, varying from gamma1 to gamma1end
    
    if Lpz(1,1)>L                   % sk=[2] exists only if lpz[2]>L
        gamma1=gmIn(1);
        gamma1end=gmIn(2)-tau1*L^2/(T*tm*Eb);
        gammaStepNumber=100;
        gammaStep=(gamma1end-gamma1)/gammaStepNumber;
        gamma=gamma1:gammaStep:gamma1end
        
        sigma2=repmat(tau1*L/T,1,gammaStepNumber+1)
        epsilon2=100*(3/2*tau1*L/(T*Eb)+gamma*tm/(2*L))
        
        hold on
        plot(epsilon2,sigma2,'g')
        
        hold on
        plot(epsilon2(end),sigma2(end),'g*')

    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Find theoritical function between epsilon and sigma sk=[1,2,3] %%%%%%
        
    if Lpz(1,1)<L                   % sk=[1,2,3] exists only if lpz[2]<L

        

    end    
    
end
end

