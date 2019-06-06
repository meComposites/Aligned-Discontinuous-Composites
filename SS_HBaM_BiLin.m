function [EInt,deInt,NewEpsilon,NewSigma] = SS_HBaM_BiLin(EE,DE,ET,DT,tm,L,SmIn,Giic,G)
%
%%%% As provided by Joel, but with corrected P1_F_123 and P1_FC_123
%%%% Modified on Sunday 13/08/2017 at 17:54: 
%%%%    Removed post-processing of the curve;
%%%%    Strains in absolute values (not in %)
    
% L     Length of the overlap
% e.g.
% SS_HBaM(933000,587000,0.017/4,0.003/4,0.01,1.5,80,1,1500)
% Path 1: SS_HBaM_BiLin(933000,787000,0.017/4,0.003/4,0.01,1.5,80,1,1500)
% Path 2: SS_HBaM_BiLin(933000,787000,0.017/4,0.003/4,0.01,0.1,80,1,1500)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define aux variables
% Search parameters
opts = optimset('TolX',10^(-32),'Display','off');

% Constitutive law
gm=[SmIn/G SmIn/G+(2*Giic/tm-(SmIn*SmIn/G))/SmIn 10^5]';
Sm=[SmIn 0 0]';
G=[G (Sm(2)-Sm(1))/(gm(2)-gm(1)) 0]';

% Parameters
for ig=1:1:length(G) % Lambda
    Y(ig)=sqrt(8*abs(G(ig))/tm*(EE*ET+DE*DT)/((EE^2-DE^2)*(ET^2-DT^2)));
end
K=2*(DE*ET)/(EE*ET+DE*DT); % Kappa
p=(2*ET/(ET-DT)+K)/(2*ET/(ET+DT)-K); % Rho (EB.TE)/(EA.TA) ??
p=max(p,1/p); % to get the stiffest on the same size?
type=sign(G); % Type of subdomains
N=length(G)-1; % Number of subdomains
Teq=((ET^2-DT^2)/(2*ET)); % Equivalent thickness
AR=0;AL=0;g_act=0;
Linear_SS=F_EndLinear(DE,DT,EE,ET,tm,L/2,G(1));


%%
% intermediate points
goi=0;
ni=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path identification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Paths:         Path 1      Path 2      Path 3
%               [1]         [1]         [1]
%               [1 2]       [1 2]       [1 2]
%               [1 2 3]     [2 1 2]     [2 1 2]
%               [2 1 2 3]   [2]         [2 1 2 3]
%               [2 3]       [3]         [2 3]
%               [3]                     [3]

I1_PsiL=2*SmIn/(Y(2)*Teq);
I1_PsiR=-I1_PsiL/p;
I1_tR=sqrt(SmIn^2-I1_PsiR^2*Y(2)^2*Teq^2/4);
I1_lt=2/Y(2)*(atan(1)+atan(sqrt((SmIn-I1_tR)/(SmIn+I1_tR))));

if I1_lt>L
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Path 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I_Path=2;
    P2_SS=Linear_SS;
%%    
    % find [2 1 2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find go and calculate SS response:
    [P2_g0_212,jkl,~]= fminbnd(@P2_F_212,0,gm(1),opts); % Good
    g_act=[g_act P2_g0_212];
    [PsiL,gL,tL] = P2_FC_212(P2_g0_212);
    [epsilon_inf, sigma_inf] = post_process(PsiL,gL,gm(1));
    % calculate ni intermediate SS point:
    % ni=1;
    goint=exp(log(g_act(end-1)):(log(g_act(end))-log(g_act(end-1)))/ni:log(g_act(end)));
    goint=goint(2:end-1);
    for noi=1:length(goint)
        goi=goint(noi);
        [P2_tL_in_12,jkl,~]= fminbnd(@P2_FC_errin_12,tL,SmIn,opts); % Good
        [PsiR2,p2_gL,p2_gR]= P2_FC_in_12(P2_tL_in_12);
        [epsilon_infi, sigma_infi] = post_process(abs(PsiR2),p2_gL,p2_gR);
        P2_SS=[P2_SS;epsilon_infi sigma_infi];
    end
    % Update SS curve
    P2_SS=[P2_SS;epsilon_inf sigma_inf];
    
    % find [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find tL and calculate SS response:
    [P2_tL_2,jkl,~]= fminbnd(@P2_F_2,0,tL,opts); % Good
    g_act=[g_act gm(1)];
    [PsiR,gL,gR,P2_lL]=P2_FC_2(P2_tL_2);
    [epsilon_inf, sigma_inf] = post_process(abs(PsiR),gL,gR);
    % calculate ni intermediate SS point:
    goint=exp(log(g_act(end-1)):(log(g_act(end))-log(g_act(end-1)))/ni:log(g_act(end)));
    goint=goint(2:end-1);
    for noi=1:length(goint)
        goi=goint(noi);
        [P2_tL_in_12,jkl,~]= fminbnd(@P2_FC_errin_212,P2_tL_2,tL,opts); % Good
        [PsiR2,p2_gL,p2_gR]= P2_FC_in_212(P2_tL_in_12);
        [epsilon_infi, sigma_infi] = post_process(abs(PsiR2),p2_gL,p2_gR);
        P2_SS=[P2_SS;epsilon_infi sigma_infi];
    end
    % Update SS curve
    
    P2_SS=[P2_SS;epsilon_inf sigma_inf];
    
    % find [3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g_act=[g_act gm(2)];
    [PsiR,gL,gR] = P2_FC_3(P2_lL);
    [epsilon_inf, sigma_inf] = post_process(abs(PsiR),gL,gR);
    % Update SS curve
    P2_SS=[P2_SS;epsilon_inf sigma_inf];
    
%     figure
%     hold on
%     plot(P2_SS(:,1),P2_SS(:,2))
        
    Epsilon=P2_SS(:,1);
    Sigma=P2_SS(:,2);
    EInt=Sigma(2)/Epsilon(2);
    deInt=Epsilon(end)-Epsilon(end-1);

else
    %[I2_eps,jkl,~]= fminsearch(@I2_F_eps,4.3762e-31,opts); % Not great
    I2_t0_min=SmIn/(cosh(Y(1)*L/2));
    
    [I2_t0G_sec,jkl,~]= fminbnd(@I2_F_eps_sec,I2_t0_min,SmIn,opts); % Good
    %I2_Fcrack_sec(I2_t0G_sec)
    
%     [I2_t0G,jkl,~]= fminbnd(@I2_F_eps,I2_t0_min,SmIn,opts); % Good
%     I2_tRf=I2_F_tR(I2_t0G);
    if ~I2_Fcrack_sec(I2_t0G_sec) % if there is no crack
%     if I2_tRf>SmIn && jkl <0.001
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Path 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I_Path=3;
        P3_SS=Linear_SS;
        
        % find [2 1 2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find go and calculate SS response:
        [P3_g0_212,jkl,~]= fminbnd(@P2_F_212,0,gm(1),opts); % Good
        g_act=[g_act P3_g0_212];
        [PsiL,gL,tL] = P2_FC_212(P3_g0_212);
        [epsilon_inf, sigma_inf] = post_process(PsiL,gL,gm(1));
        % calculate ni intermediate SS point:
        goint=exp(log(g_act(end-1)):(log(g_act(end))-log(g_act(end-1)))/ni:log(g_act(end)));
        goint=goint(2:end-1);
        for noi=1:length(goint)
            goi=goint(noi);
            [P3_tL_in_12,jkl,~]= fminbnd(@P2_FC_errin_12,tL,SmIn,opts); % Good
            [PsiR2,p2_gL,p2_gR]= P2_FC_in_12(P3_tL_in_12);
            [epsilon_infi, sigma_infi] = post_process(abs(PsiR2),p2_gL,p2_gR);
            P3_SS=[P3_SS;epsilon_infi sigma_infi];
        end
        % Update SS curve
        P3_SS=[P3_SS;epsilon_inf sigma_inf];
        
        % find [2 1 2 3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find go and calculate SS response:
        [P3_g0_2123,jkl,~]= fminbnd(@P3_F_2123,g_act(end),gm(1),opts); % Good
        g_act=[g_act P3_g0_2123];
        [PsiL,gL,gR] = P2_FC_2123(P3_g0_2123);
        [epsilon_inf, sigma_inf] = post_process(abs(PsiL),gL,gR);
        % calculate ni intermediate SS point:
        goint=exp(log(g_act(end-1)):(log(g_act(end))-log(g_act(end-1)))/ni:log(g_act(end)));
        goint=goint(2:end-1);
        for noi=1:length(goint)
            goi=goint(noi);
            [P3_tL_in_212,jkl,~]= fminbnd(@P2_FC_errin_212,0,tL,opts); % Good
            % [P3_tL_in_212,jkl,~]= fminbnd(@P2_FC_errin_12,tL,SmIn,opts); % Good
            [PsiR2,p2_gL,p2_gR]= P2_FC_in_12(P3_tL_in_212);
            [epsilon_infi, sigma_infi] = post_process(abs(PsiR2),p2_gL,p2_gR);
            P3_SS=[P3_SS;epsilon_infi sigma_infi];
        end
        P3_SS=[P3_SS;epsilon_inf sigma_inf];
        
        % find [2 3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        g_act=[g_act gm(2)];
        [PsiL,gL,gR] = P3_FC_2123();
        [epsilon_inf, sigma_inf] = post_process(abs(PsiL),gL,gR);
        P3_SS=[P3_SS;epsilon_inf sigma_inf];
        

        
        % find [3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         figure
%         hold on
%         plot(P3_SS(:,1),P3_SS(:,2))
        
        Epsilon=P3_SS(:,1);
        Sigma=P3_SS(:,2);
        EInt=Sigma(2)/Epsilon(2);
        deInt=0;
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Path 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I_Path=1;
        P1_SS=Linear_SS;
        
        % find [1 2 3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find go and calculate SS response:
        [P1_g0_123,jkl,~]= fminbnd(@P1_F_123,g_act(end),gm(1),opts); % Good
        g_act=[g_act P1_g0_123];
        [PsiL,gL,gR] = P1_FC_123(P1_g0_123);
        [epsilon_inf, sigma_inf] = post_process(abs(PsiL),gL,gR);
        % calculate ni intermediate SS point:
        goint=exp(log(g_act(end-1)):(log(g_act(end))-log(g_act(end-1)))/ni:log(g_act(end)));
        goint=goint(2:end-1);
        for noi=1:length(goint)
            goi=goint(noi);
            [P1_tL_in_212,jkl,~]= fminbnd(@P1_FC_errin_212,0,SmIn,opts); % Good
            % [P3_tL_in_212,jkl,~]= fminbnd(@P2_FC_errin_12,tL,SmIn,opts); % Good
            [PsiR2,p1_gL,p1_gR]= P1_FC_in_12(P1_tL_in_212);
            [epsilon_infi, sigma_infi] = post_process(abs(PsiR2),p1_gL,p1_gR);
            P1_SS=[P1_SS;epsilon_infi sigma_infi];
        end
        P1_SS=[P1_SS;epsilon_inf sigma_inf];
        
        % find [2 1 2 3] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find go and calculate SS response:
        [P1_g0_2123,jkl,~]= fminbnd(@P1_F_2123,g_act(end),gm(1),opts); % Good
        g_act=[g_act P1_g0_2123];
        [PsiL,gL,gR] = P1_FC_2123(P1_g0_2123);
        [epsilon_inf, sigma_inf] = post_process(abs(PsiL),gL,gR);
        % calculate 1 intermediate SS point: (No need because flat?)
%         [P3_tL_in_123,jkl,~]= fminbnd(@P2_FC_errin_212,0,tL,opts); % Good
%         [PsiR2,p2_gL,p2_gR]= P2_FC_in_12(P3_tL_in_123);
%         [epsilon_inf2, sigma_inf2] = post_process(abs(PsiR2),p2_gL,p2_gR);
        P1_SS=[P1_SS;epsilon_inf sigma_inf];
        
%         figure
%         hold on
%         plot(P1_SS(:,1),P1_SS(:,2))
        
        Epsilon=P1_SS(:,1);
        Sigma=P1_SS(:,2);
        EInt=Sigma(2)/Epsilon(2);
        deInt=0;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Npercents=1.5;  % (%)

% Correction for plateau
[M,indBB]=max(Sigma); % Sigma Max
indMin=min(find(Sigma>max(Sigma)*(1-Npercents/100)));
emin=Epsilon(indMin); % epsilon min in N%
smin=max(Sigma)*(1-Npercents/100);
smin=Sigma(indMin);
indMax=max(find(Sigma>max(Sigma)*(1-Npercents/100)));
emax=max(max(Epsilon(indMin:indMax)),Epsilon(indMax)); % epsilon max in N%
smax=max(Sigma);

% New curve with loading only
if smin==smax
    NewSigma=[Sigma(1:indMin-1) ; smax];
    NewEpsilon=[Epsilon(1:indMin-1) ; emax];
else
    NewSigma=[Sigma(1:indMin-1) ; smin ; smax];
    NewEpsilon=[Epsilon(1:indMin-1) ; emin ; emax];
end


return;

    function [epsilon_inf, sigma_inf] = post_process(PsiR,gL,gR)
        if AL==0
            sigma_inf=PsiR/(2*ET/(ET-DT)-K);
            epsilon_inf=(tm/(2*L)*((gL*(EE+DE)*(ET-DT)+gR*(EE-DE)*(ET+DT))/(EE*ET+DE*DT))+2*sigma_inf*ET/(EE*ET+DE*DT));
        else
            sigma_inf=PsiR/(2*ET/(ET+DT)-K);
            epsilon_inf=(tm/(2*L)*((gL*(EE+DE)*(ET+DT)+gR*(EE-DE)*(ET-DT))/(EE*ET+DE*DT))+2*sigma_inf*ET/(EE*ET+DE*DT));
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions identification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Linear_SS = F_EndLinear(DE,DT,EE,ET,tm,l,G)
        [~,Lchar,sig_i,eps_i] = Eini_2Stiffness(DE,DT,EE,ET,tm,l,G);
        if L-Lchar > Lchar       % The right side is the longest
            g_act = gm(1)/cosh(Y(1)*(L-Lchar));
            AR=1;
            AL=0;
        else                       % The left side is the longest
            g_act = gm(1)/cosh(Y(1)*Lchar);
            AR=0;
            AL=1;
        end
        g_i=[0;g_act];
        sig_ini=sig_i*G(1).*g_i;
        eps_ini=eps_i*G(1).*g_i;
        Linear_SS=[eps_ini sig_ini];
    end
    function I2_eps = I2_F_eps(t0_guess)
        I2_PsiL1=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0_guess^2);
        I2_l1=2/Y(1)*atanh(sqrt((SmIn-t0_guess)/(SmIn+t0_guess)));
        I2_l1(I2_l1==Inf)=1/Y(1)*log(2*SmIn/t0_guess);
        
        I2_PsiL=-sqrt(I2_PsiL1^2+(2/(Y(2)*Teq))^2*SmIn^2);
        I2_l2=2/Y(2)*atan(1/SmIn*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0_guess^2)+SmIn^2)+(Y(2)*Teq/2/SmIn)*I2_PsiL1);
        
        I2_th2=tanh(Y(1)/2*(L-I2_l1-I2_l2))^2;
        I2_tR=t0_guess*(1+I2_th2)/(1-I2_th2);
        I2_tR(I2_tR==Inf)=2*t0_guess*exp(Y(1)*(L-I2_l1-I2_l2));
        
        I2_tR=t0_guess*cosh(Y(1)*(L-I2_l1-I2_l2));
        
        I2_PsiR=2/(Y(1)*Teq)*sqrt(I2_tR^2-t0_guess^2);
        
        I2_eps=abs(I2_PsiL+1/p*I2_PsiR);
    end

    function I2_eps = I2_F_eps_sec(t0_guess)
        I2_PsiL1=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0_guess^2);
        I2_l1=2/Y(1)*atanh(sqrt((SmIn-t0_guess)/(SmIn+t0_guess)));
        I2_l1(I2_l1==Inf)=1/Y(1)*log(2*SmIn/t0_guess);
        
        I2_l2=2/Y(2)*atan(1/SmIn*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0_guess^2)+SmIn^2)+(Y(2)*Teq/2/SmIn)*I2_PsiL1);
        
        if I2_l2>L-2*I2_l1 % if can't reach the crack
            I2_l2=L-2*I2_l1;
            I2_tL=-I2_PsiL1*Y(2)*Teq/2*sin(Y(2)*I2_l2)+SmIn*cos((Y(2)*I2_l2));
            I2_PsiL=-sqrt(I2_PsiL1^2+(2/(Y(2)*Teq))^2*(SmIn^2-I2_tL^2));
        else % if can reach the crack
            I2_PsiL=-sqrt(I2_PsiL1^2+(2/(Y(2)*Teq))^2*SmIn^2);
        end
        I2_PsiR=-I2_PsiL1;
        
        I2_eps=abs(I2_PsiL+p*I2_PsiR);
    end

    function I2_crack = I2_Fcrack_sec(t0)
        I2_PsiL1=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0^2);
        I2_l1=2/Y(1)*atanh(sqrt((SmIn-t0)/(SmIn+t0)));
        I2_l1(I2_l1==Inf)=1/Y(1)*log(2*SmIn/t0);
        
        I2_l2=2/Y(2)*atan(1/SmIn*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0^2)+SmIn^2)+(Y(2)*Teq/2/SmIn)*I2_PsiL1);
        
        I2_crack=L>I2_l2+2*I2_l1;
    end
    function I2_tR = I2_F_tR(t0_guess)
        I2_PsiL1=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0_guess^2);
        I2_l1=2/Y(1)*atanh(sqrt((SmIn-t0_guess)/(SmIn+t0_guess)));
        I2_l1(I2_l1==Inf)=1/Y(1)*log(2*SmIn/t0_guess);
        
        I2_PsiL=-sqrt(I2_PsiL1^2+(2/(Y(2)*Teq))^2*SmIn^2);
        I2_l2=2/Y(2)*atan(1/SmIn*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0_guess^2)+SmIn^2)+(Y(2)*Teq/2/SmIn)*I2_PsiL1);
        
        I2_th2=tanh(Y(1)/2*(L-I2_l1-I2_l2))^2;
        I2_tR=t0_guess*(1+I2_th2)/(1-I2_th2);
        I2_tR(I2_tR==Inf)=2*t0_guess*exp(Y(1)*(L-I2_l1-I2_l2));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions Path 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function P2_212 = P2_F_212(g0_guess)
        t0_guess=g0_guess*G(1);
        P2_PsiL=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0_guess^2);
        P2_l1=2/Y(1)*atanh(sqrt((SmIn-t0_guess)/(SmIn+t0_guess)));
        P2_l1(P2_l1==Inf)=1/Y(1)*log(2*SmIn/t0_guess);
        
        P2_lR=L-2*P2_l1;
        
        P2_tR=sqrt(SmIn^2-(SmIn^2-t0_guess^2)*(Y(2)/Y(1))^2*(p^2-1));
        
        P2_lR2=2/Y(2)*atan(1/(SmIn+P2_tR)*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0_guess^2)+(SmIn^2-P2_tR^2))+(Y(2)*Teq/2/(SmIn+P2_tR))*P2_PsiL);
        
        P2_212=abs(P2_lR-P2_lR2);
    end
    function [P2_PsiL,gL,P2_tL] = P2_FC_212(g0_guess)
        t0_guess=g0_guess*G(1);
        P2_PsiL=2/(Y(1)*Teq)*sqrt(SmIn^2-t0_guess^2);
        P2_tL=sqrt(SmIn^2-(SmIn^2-t0_guess^2)*(Y(2)/Y(1))^2*(p^2-1));
        gL=gm(2)-(gm(2)-gm(1))*P2_tL/SmIn;
        %PsiR=-p*P2_PsiL;
    end
    function [P2_err] = P2_FC_errin_12(tL_guess)
        g0=goi;
        t0=g0*G(1);

        P2_PsiL=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0^2);
        P2_l1=2/Y(1)*atanh(sqrt((SmIn-t0)/(SmIn+t0)));
        P2_l1(P2_l1==Inf)=1/Y(1)*log(2*SmIn/t0);
        
        P2_lL2=2/Y(2)*atan(1/(SmIn+tL_guess)*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0^2)+(SmIn^2-tL_guess^2))+(Y(2)*Teq/2/(SmIn+tL_guess))*P2_PsiL);
        P2_PsiL2=sqrt(P2_PsiL^2+(2/(Y(2)*Teq))^2*(SmIn^2-tL_guess^2));
        P2_PsiR=-P2_PsiL2/p;
        
        P2_tR=sqrt(t0^2+(Y(1)*Teq*P2_PsiR/2)^2);
        P2_lR1=2/Y(1)*atanh(sqrt((P2_tR-t0)/(P2_tR+t0)));
        P2_lR1(P2_lR1==Inf)=1/Y(1)*log(2*P2_tR/t0);
        
        P2_err=abs(L-P2_l1-P2_lL2-P2_lR1);
    end
    function [P2_PsiR,p2_gL,p2_gR] = P2_FC_in_12(tL)
        g0=goi;
        t0=g0*G(1);

        P2_PsiL=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0^2);
        P2_l1=2/Y(1)*atanh(sqrt((SmIn-t0)/(SmIn+t0)));
        P2_l1(P2_l1==Inf)=1/Y(1)*log(2*SmIn/t0);
        
        P2_lL2=2/Y(2)*atan(1/(SmIn+tL)*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0^2)+(SmIn^2-tL^2))+(Y(2)*Teq/2/(SmIn+tL))*P2_PsiL);
        P2_PsiL2=sqrt(P2_PsiL^2+(2/(Y(2)*Teq))^2*(SmIn^2-tL^2));
        P2_PsiR=-P2_PsiL2/p;
        
        P2_tR=sqrt(t0^2+(Y(1)*Teq*P2_PsiR/2)^2);
        p2_gR=P2_tR/G(1);
        p2_gL=gm(2)-(gm(2)-gm(1))*tL/SmIn;
    end

    function [P2_err] = P2_F_2(tL_guess)

        g0=gm(1);
        t0=SmIn;

        P2_PsiL=-2/(Y(2)*Teq)*sqrt(t0^2-tL_guess^2);
        P2_l1L=2/Y(2)*atan(1*sqrt((SmIn-tL_guess)/(SmIn+tL_guess)));
        
        P2_l1R=L-P2_l1L;
        P2_PsiR=-P2_PsiL/p;
        P2_tR=sqrt(t0^2-(Y(2)*Teq*P2_PsiR/2)^2);
        P2_lR1=2/Y(2)*atanh(sqrt((t0-P2_tR)/(t0+P2_tR)));
        
        P2_err=abs(P2_l1R-P2_lR1);
    end
    function [P2_PsiR,gL,gR,P2_l1L] = P2_FC_2(tL)
        g0=gm(1);
        t0=SmIn;
        gL=gm(2)-(gm(2)-gm(1))*tL/SmIn;
        
        P2_PsiL=-2/(Y(2)*Teq)*sqrt(t0^2-tL^2);
        P2_l1L=2/Y(2)*atan(1*sqrt((SmIn-tL)/(SmIn+tL)));
        
        P2_PsiR=-P2_PsiL/p;
        P2_tR=sqrt(t0^2-(Y(2)*Teq*P2_PsiR/2)^2);
        gR=gm(2)-(gm(2)-gm(1))*P2_tR/SmIn;
    end
    function [P2_err] = P2_FC_errin_212(tL_guess)
        g0=goi;
        t0=g0*G(1);
        
        P2_PsiL=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0^2);
        P2_l1=2/Y(1)*atanh(sqrt((SmIn-t0)/(SmIn+t0)));
        P2_l1(P2_l1==Inf)=1/Y(1)*log(2*SmIn/t0);
        
        P2_lL2=2/Y(2)*atan(1/(SmIn+tL_guess)*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0^2)+(SmIn^2-tL_guess^2))+(Y(2)*Teq/2/(SmIn+tL_guess))*P2_PsiL);
        
        P2_lR2=L-2*P2_l1-P2_lL2;
        
        P2_PsiL2=sqrt(P2_PsiL^2+(2/(Y(2)*Teq))^2*(SmIn^2-tL_guess^2));
        P2_PsiR2=-P2_PsiL2/p;
        
        P2_tR=sqrt(SmIn^2-(P2_PsiR2^2-P2_PsiL^2)*((Y(2)*Teq)/2)^2);
        P2_lR22=2/Y(2)*atan(1/(SmIn+P2_tR)*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0^2)+(SmIn^2-P2_tR^2))+(Y(2)*Teq/2/(SmIn+P2_tR))*P2_PsiL);
        
        P2_err=abs(P2_lR2-P2_lR22);

    end
    function [P2_PsiR2,p2_gL,p2_gR] = P2_FC_in_212(tL)
        g0=goi;
        t0=g0*G(1);
        
        P2_PsiL=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0^2);
        P2_PsiL2=sqrt(P2_PsiL^2+(2/(Y(2)*Teq))^2*(SmIn^2-tL^2));
        P2_PsiR2=-P2_PsiL2/p;
        
        P2_tR=sqrt(SmIn^2-(P2_PsiR2^2-P2_PsiL^2)*((Y(2)*Teq)/2)^2);
        
        p2_gL=gm(2)-(gm(2)-gm(1))*tL/SmIn;
        p2_gR=gm(2)-(gm(2)-gm(1))*P2_tR/SmIn;
    end

    function [P2_PsiR,gL,gR] = P2_FC_3(P2_lL)
        g0=gm(2)*.999;
        t0=SmIn*(1-(g0-gm(1))/(gm(2)-gm(1)));
        
        tL=t0*cos(Y(2)*P2_lL);
        gL=gm(2)-(gm(2)-gm(1))*tL/SmIn;
        
        P2_PsiL=-2/(Y(2)*Teq)*sqrt(t0^2-tL^2);
        P2_PsiR=-P2_PsiL/p;
        P2_tR=sqrt(t0^2-(Y(2)*Teq*P2_PsiR/2)^2);
        gR=gm(2)-(gm(2)-gm(1))*P2_tR/SmIn;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions Path 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function P_err_123 = P1_F_123(g0_guess)
t0_guess=g0_guess*G(1);
        P1_PsiL=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0_guess^2);
        P1_l1=2/Y(1)*atanh(sqrt((SmIn-t0_guess)/(SmIn+t0_guess)));
        P1_l1(P1_l1==Inf)=1/Y(1)*log(2*SmIn/t0_guess);
        
        P1_lL2=2/Y(2)*atan(1/(SmIn)*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0_guess^2)+(SmIn^2))+(Y(2)*Teq/2/(SmIn))*P1_PsiL);
        P1_PsiL2=sqrt(P1_PsiL^2+(2/(Y(2)*Teq))^2*(SmIn^2));
        
        P1_PsiR2=-P1_PsiL2/p;
        
        %P1_tR=sqrt(SmIn^2-((Y(1)*Teq*P1_PsiR2)/2)^2);
        P1_tR=sqrt(t0_guess^2+((Y(1)*Teq*P1_PsiR2)/2)^2);
        P1_lR2=2/Y(1)*atanh(sqrt((P1_tR-t0_guess)/(P1_tR+t0_guess)));
        
        P_err_123=abs(L-P1_l1-P1_lL2-P1_lR2);
    end
    function [P1_PsiR2,P1_gL,P1_gR] = P1_FC_123(g0)
        t0=g0*G(1);
        P1_PsiL=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0^2);
        P1_l1=2/Y(1)*atanh(sqrt((SmIn-t0)/(SmIn+t0)));
        P1_l1(P1_l1==Inf)=1/Y(1)*log(2*SmIn/t0);
        
        P1_PsiL2=sqrt(P1_PsiL^2+(2/(Y(2)*Teq))^2*(SmIn^2));
        P1_PsiR2=-P1_PsiL2/p;
        
        %P1_tR=sqrt(SmIn^2-((Y(1)*Teq*P1_PsiR2)/2)^2);
        P1_tR=sqrt(t0^2+((Y(1)*Teq*P1_PsiR2)/2)^2);
        P1_gR=gm(1)*P1_tR/SmIn;
        P1_gL=gm(2);
       
    end
    function P1_err_2123 = P1_F_2123(g0_guess)
        t0_guess=g0_guess*G(1);
        P1_PsiL=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0_guess^2);
        P1_l1=2/Y(1)*atanh(sqrt((SmIn-t0_guess)/(SmIn+t0_guess)));
        P1_l1(P1_l1==Inf)=1/Y(1)*log(2*SmIn/t0_guess);
        
        P1_lL2=2/Y(2)*atan(1/(SmIn)*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0_guess^2)+(SmIn^2))+(Y(2)*Teq/2/(SmIn))*P1_PsiL);
        P1_PsiL2=sqrt(P1_PsiL^2+(2/(Y(2)*Teq))^2*(SmIn^2));
        
        P1_PsiR2=-P1_PsiL2/p;
        
%         
%         P3_tR=sqrt(SmIn^2-(P3_PsiR2^2-P3_PsiL^2)*((Y(2)*Teq)/2)^2);
%         P3_lR2=2/Y(2)*atan(1/(SmIn+P3_tR)*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0_guess^2)+(SmIn^2-P3_tR^2))+(Y(2)*Teq/2/(SmIn+P3_tR))*P3_PsiL);
        
        P1_err_2123=abs(P1_PsiR2-P1_PsiL);
    end
    function [P1_PsiR2,P1_gL,P1_gR] = P1_FC_2123(g0_guess)
        t0_guess=g0_guess*G(1);
        P1_PsiL=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0_guess^2);
        P1_l1=2/Y(1)*atanh(sqrt((SmIn-t0_guess)/(SmIn+t0_guess)));
        P1_l1(P1_l1==Inf)=1/Y(1)*log(2*SmIn/t0_guess);
        
        P1_lL2=2/Y(2)*atan(1/(SmIn)*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0_guess^2)+(SmIn^2))+(Y(2)*Teq/2/(SmIn))*P1_PsiL);
        P1_PsiL2=sqrt(P1_PsiL^2+(2/(Y(2)*Teq))^2*(SmIn^2));
        
        P1_PsiR2=-P1_PsiL2/p;
        
        P1_Ldeb=L-2*P1_l1-P1_lL2;
        P1_gL=gm(2)+P1_PsiL2*P1_Ldeb*2/(tm*ET)*(EE*ET+DE*DT)/(EE^2-DE^2);
        P1_gR=gm(1);
    end

    function [P2_PsiR,p2_gL,p2_gR] = P1_FC_in_12(tL)
        g0=goi;
        t0=g0*G(1);

        P2_PsiL=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0^2);
        P2_l1=2/Y(1)*atanh(sqrt((SmIn-t0)/(SmIn+t0)));
        P2_l1(P2_l1==Inf)=1/Y(1)*log(2*SmIn/t0);
        
        P2_lL2=2/Y(2)*atan(1/(SmIn+tL)*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0^2)+(SmIn^2-tL^2))+(Y(2)*Teq/2/(SmIn+tL))*P2_PsiL);
        P2_PsiL2=sqrt(P2_PsiL^2+(2/(Y(2)*Teq))^2*(SmIn^2-tL^2));
        P2_PsiR=-P2_PsiL2/p;
        
        P2_tR=sqrt(t0^2+(Y(1)*Teq*P2_PsiR/2)^2);
        p2_gR=P2_tR/G(1);
        p2_gL=gm(2)-(gm(2)-gm(1))*tL/SmIn;
    end
    function [P2_err] = P1_FC_errin_212(tL_guess)
        g0=goi;
        t0=g0*G(1);
        
        P2_PsiL=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0^2);
        P2_l1=2/Y(1)*atanh(sqrt((SmIn-t0)/(SmIn+t0)));
        P2_l1(P2_l1==Inf)=1/Y(1)*log(2*SmIn/t0);
        
        P2_lL2=2/Y(2)*atan(1/(SmIn+tL_guess)*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0^2)+(SmIn^2-tL_guess^2))+(Y(2)*Teq/2/(SmIn+tL_guess))*P2_PsiL);
        
        P2_lR2=L-2*P2_l1-P2_lL2;
        
        P2_PsiL2=sqrt(P2_PsiL^2+(2/(Y(2)*Teq))^2*(SmIn^2-tL_guess^2));
        P2_PsiR2=-P2_PsiL2/p;
        
        P2_tR=sqrt(SmIn^2-(P2_PsiR2^2-P2_PsiL^2)*((Y(2)*Teq)/2)^2);
        P2_lR22=2/Y(2)*atan(1/(SmIn+P2_tR)*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0^2)+(SmIn^2-P2_tR^2))+(Y(2)*Teq/2/(SmIn+P2_tR))*P2_PsiL);
        
        P2_err=abs(P2_lR2-P2_lR22);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions Path 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function P3_err_2123 = P3_F_2123(g0_guess)
        t0_guess=g0_guess*G(1);
        P3_PsiL=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0_guess^2);
        P3_l1=2/Y(1)*atanh(sqrt((SmIn-t0_guess)/(SmIn+t0_guess)));
        P3_l1(P3_l1==Inf)=1/Y(1)*log(2*SmIn/t0_guess);
        
        P3_lL2=2/Y(2)*atan(1/(SmIn)*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0_guess^2)+(SmIn^2))+(Y(2)*Teq/2/(SmIn))*P3_PsiL);
        P3_PsiL2=sqrt(P3_PsiL^2+(2/(Y(2)*Teq))^2*(SmIn^2));
        
        P3_PsiR2=-P3_PsiL2/p;
        
        P3_tR=sqrt(SmIn^2-(P3_PsiR2^2-P3_PsiL^2)*((Y(2)*Teq)/2)^2);
        P3_lR2=2/Y(2)*atan(1/(SmIn+P3_tR)*sqrt((Y(2)/Y(1))^2*(SmIn^2-t0_guess^2)+(SmIn^2-P3_tR^2))+(Y(2)*Teq/2/(SmIn+P3_tR))*P3_PsiL);
        
        P3_err_2123=abs(L-2*P3_l1-P3_lL2-P3_lR2);
    end
    function [P3_PsiR2,P3_gL,P3_gR] = P2_FC_2123(g0)
        t0=g0*G(1);
        P3_PsiL=-2/(Y(1)*Teq)*sqrt(SmIn^2-t0^2);
        P3_l1=2/Y(1)*atanh(sqrt((SmIn-t0)/(SmIn+t0)));
        P3_l1(P3_l1==Inf)=1/Y(1)*log(2*SmIn/t0);
        
        P3_PsiL2=sqrt(P3_PsiL^2+(2/(Y(2)*Teq))^2*(SmIn^2));
        P3_PsiR2=-P3_PsiL2/p;
        P3_tR=sqrt(SmIn^2-(P3_PsiR2^2-P3_PsiL^2)*((Y(2)*Teq)/2)^2);
        P3_gR=gm(2)-(gm(2)-gm(1))*P3_tR/SmIn;
        P3_gL=gm(2);
       
    end
    function [P3_C_PsiR,gL,gR] = P3_FC_2123()
        P3_C_PsiL=2*SmIn/(Y(2)*Teq);
        P3_C_PsiR=-P3_C_PsiL/p;
        P3_C_tR=sqrt(SmIn^2-P3_C_PsiR^2*Y(2)^2*Teq^2/4);
        
        gR=gm(2)-(gm(2)-gm(1))*P3_C_tR/SmIn;
        I1_C_lt=2/Y(2)*(atan(1)+atan(sqrt((SmIn-P3_C_tR)/(SmIn+P3_C_tR))));
        P3_C_Ldeb=L-I1_C_lt;
        
        gL=gm(2)+P3_C_PsiL*P3_C_Ldeb*2/(tm*ET)*(EE*ET+DE*DT)/(EE^2-DE^2);
    end



end

