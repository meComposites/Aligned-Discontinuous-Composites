function [ SS_curve,EEE,strength,eps_softening,SS_steps_R,Giic,vec_l_int ] = general_loop_2Stiffness_nonIterative_VariableGiicSm(EE,DE,ET,DT,tm,L,Giic,MLawS,MLawg)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts = optimset('TolX',10^(-32),'Display','off');

ng=20;
nx=10;

S3Dplot=[];
P3Dplot=[];
Xa3Dplot=[];
Xb3Dplot=[];

[G,Sm,gm,Y] = inputs_2Stiffness_VariableGiicSm(EE,DE,ET,DT,tm,Giic,MLawS,MLawg);

K=2*(DE*ET)/(EE*ET+DE*DT);
p=(2*ET/(ET-DT)+K)/(2*ET/(ET+DT)-K);

% Definition of Matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type=sign(G);
N=length(G)-1;
[Lpz,P0pz] = LpzMatrix(N,type,Y,gm,Sm,ET,DT,tm,EE,DE);
[gpzMax] = gMatrix(N,type,Y,gm,Sm,ET,DT,tm,EE,DE);
LnR=-1;

global l_int l_char_ini Xinf_initial_guess;
global DX_L X_L g_L S_L Dl_L l_L Dli_L sol_is_zero_L;
global DX_R X_R g_R S_R Dl_R l_R Dli_R sol_is_zero_R;
global gL_end SL_end gR_end SR_end;
global lengthR lengthL;
global LRA LRB;

global State_R_previous_previous State_L_previous_previous;
global State_R_previous State_L_previous State_R State_L;

State_R_previous=1;
State_L_previous=1;

global P0iLPlot P0iRPlot g0iLPlot g0iRPlot S0iLPlot S0iRPlot lL0Plot lR0Plot;

global SS_steps;
SS_steps=[0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%                           %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%   LINEAR-ELASTIC DOMAIN   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%                           %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[EEE,Lchar,sig_i,eps_i] = Eini_2Stiffness(DE,DT,EE,ET,tm,L/4,G(1));
EEE=EEE*(ET+tm)/ET;

if L/2-Lchar > Lchar       % The right side is the longest
    g_act = gm(1)/cosh(Y(1)*(L/2-Lchar));
    State_R=[1 2];
    State_L=1;
else                       % The left side is the longest
    g_act = gm(1)/cosh(Y(1)*Lchar);
    State_L=[1 2];
    State_R=1;
end

g_i=(0:g_act/5:g_act)';
sig_ini=sig_i*G(1).*g_i;
eps_ini=eps_i*G(1).*g_i*100;

vec_l_int=[Lchar Lchar Lchar Lchar];
SS_curve=[eps_ini sig_ini];
SS_steps=SS_curve;
g_crit=[g_act];
LRB=L/2-Lchar;
States{1}=[1];
% States{2}=[fliplr(State_L) State_R(2:end)];

%     figure(8)
%     hold on;
%     plot(eps_ini(end),sig_ini(end),'*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   REST OF THE CURVE   %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while ~((State_L(end)==N+1 && State_R(end)==N+1) || ...                         %[ 3 | 3 ]
        (State_L(1)==N && State_L(end)==N+1 && State_R(end)==N+1) || ...        %[ 3 | 2 | 3 ]
        (State_L(end)==N+1 && State_R(1)==N && State_R(end)==N+1))              %[ 3 | 2 | 3 ]
    
    States{end+1}=[fliplr(State_L) State_R(2:end)];
    States{end};
    
    i=State_L(1);
    nL=State_L(end)-i;
    nR=State_R(end)-i;
    
    State_R_previous_previous=State_R_previous;
    State_L_previous_previous=State_L_previous;
    
    State_R_previous=State_R;
    State_L_previous=State_L;
    
    LRA=LRB;
    
    if nL == 0 || nR == 0  % Only one subdomain on one side
        % No possible deactivation of the center subdomain
        
        ginf=g_crit(end);
        gsup=gpzMax(i,i+nR,i+nL);

        find_Psi_end(ginf)
        RightActivationInf = (sign(PsiL) == sign(PsiL+p*PsiR));
        find_Psi_end(gsup)
        RightActivationSup = (sign(PsiL) == sign(PsiL+p*PsiR));
        
        if RightActivationInf == RightActivationSup
        RightActivation = (sign(PsiL) == sign(PsiL+p*PsiR));
        LeftActivation = 1-RightActivation;
        else
            1;
            %%% DO SOMETHING !!! %%%
            RightActivation = (sign(PsiL) == sign(PsiL+p*PsiR));
            LeftActivation = 1-RightActivation;
        end
        
        
        
        
        if (State_L(end)==N) && Lpz(min(end,State_L(1)),max(1,end+1-State_L(1))) > L/2 % if length of last subd > length of unit-cell, no activation of N+1 possible
            RightActivation=1;
            LeftActivation = 1-RightActivation;
        elseif (State_R(end)==N) && Lpz(min(end,State_R(1)),max(1,end+1-State_R(1))) > L/2 % if length of last subd > length of unit-cell, no activation of N+1 possible
            RightActivation=0;
            LeftActivation = 1-RightActivation;
        end
        
        if (State_L(end)==N+1)
            RightActivation=1;
            LeftActivation = 1-RightActivation;
        elseif (State_R(end)==N+1)
            RightActivation=0;
            LeftActivation = 1-RightActivation;
        end
   
        if LeftActivation
            State_L=[State_L (State_L(end)+1)];
            
            g0=[];
            gsup=gpzMax(i,State_R(end),State_L(end));
            if (size(gpzMax,2)~=max(State_L(end),State_R(end)))
                [g01,jkl1,~]= fminsearch(@find_g0_act_L,ginf,opts);
                [g0,jkl,~]=fminbnd(@find_g0_act_L,ginf,gsup,opts);
                if g01 > ginf && g01 < gsup && abs(g01-g0)>abs(g01)/20 && jkl<5
                    [g0,jkl,~]= fminsearch(@find_g0_act_L,ginf,opts);
                end
            else
                jkl=10;
            end
            
                    if isempty(g0) || jkl>5
                        [g0,jkl,~]= fminsearch(@find_g0_act_L,ginf,opts);
                        if g0 < ginf || jkl>5
                            divGsup=1;
                            jkl=1;
                            while jkl>0.1
                                divGsup=divGsup*2;
                                [g0,jkl,~]= fminbnd(@find_g0_act_L,ginf,ginf*divGsup,opts);
                            end
                        end
                    end
                    
            gR=gR_end;
            gL=gm(State_L(end)-1);
        elseif RightActivation
            State_R=[State_R (State_R(end)+1)];
            
            
            g0=[];
            gsup=gpzMax(i,State_R(end),State_L(end));
            if (size(gpzMax,2)~=max(State_L(end),State_R(end)))
                [g01,jkl1,~]= fminsearch(@find_g0_act_R,ginf,opts);
                [g0,jkl,~]=fminbnd(@find_g0_act_R,ginf,gsup,opts);
                if g01 > ginf && g01 < gsup && abs(g01-g0)>abs(g01)/20 && jkl<5
                    [g0,jkl,~]= fminsearch(@find_g0_act_R,ginf,opts);
                end
            else
                jkl=10;
            end
            
                    if isempty(g0) || jkl>5
                        [g0,jkl,~]= fminsearch(@find_g0_act_R,ginf,opts);
                        if g0 < ginf || jkl>5
                            divGsup=1;
                            jkl=1;
                            while jkl>0.1
                                divGsup=divGsup*2;
                                [g0,jkl,~]= fminbnd(@find_g0_act_R,ginf,ginf*divGsup,opts);
                            end
                        end
                    end
                    
                    
                    
            gL=gL_end;
            gR=gm(State_R(end)-1);
        end
        
        [epsilon_inf, sigma_inf] = post_process(PsiR,gL,gR);
        SS_curve=[SS_curve;epsilon_inf,sigma_inf];
        g_crit=[g_crit g0];
        
    else % Possible deactivation of the center subdomain
        
        g0=gm(i); % increase initial shear strain to reach the end of the
        % center subdomain (i just disappeared)
        
        
        if (State_L(end)==N+1) && (State_R(end)==N+1) % No possible activations: deactivation of center subdomain
            
            State_L = State_L(2:end);
            State_R = State_R(2:end);
            
            if nL==1
                Psi0L=0;
                lL0=0;
            else
                Psi0L=-P0pz(i,nL-1);
                lL0=Lpz(i,nL-1);
            end
            
            if nR==1
                Psi0R=0;
                lR0=0;
            else
                Psi0R=P0pz(i,nR-1);
                lR0=Lpz(i,nR-1);
            end
            
            [LnR,jkl,~]= fminsearch(@find_P,(L/2-(lR0+lL0))/2,opts);
            [LnR,jkl,~]= fminbnd(@find_P,0,(L/2-(lR0+lL0)),opts);
            gL=gL_end;
            gR=gR_end;
            PsiR=PR_end;
        else % Possible activation
            
            PsiL=-P0pz(min(i,end),max(1,min(end+1-i,nL)));
            %lL=Lpz(i,nL);
            
            ginf=g_crit(end);
            
            if nR==1
                PsiR1=0;
            else
                PsiR1=P0pz(i,nR-1);
            end
            PsiR2=P0pz(i,min(nR,end+1-i));
            
            %         LeftActivation = (sign(PsiL+p*PsiR1) == sign(PsiL+p*PsiR2));
            %         RightActivation = 1-LeftActivation;
            
            RightActivation = (sign(PsiL+p*PsiR1) == sign(PsiL+p*PsiR2));
            LeftActivation = 1-RightActivation;
            
        find_Psi_end(ginf)
        RightActivationInf = (sign(PsiL+p*PsiR1) == sign(PsiL+p*PsiR2));
        find_Psi_end(gsup)
        RightActivationSup = (sign(PsiL+p*PsiR1) == sign(PsiL+p*PsiR2));
        
        if RightActivationInf == RightActivationSup
        RightActivation = (sign(PsiL+p*PsiR1) == sign(PsiL+p*PsiR2));
        LeftActivation = 1-RightActivation;
        else
            1;
            %%% DO SOMETHING !!! %%%
            RightActivation = (sign(PsiL+p*PsiR1) == sign(PsiL+p*PsiR2));
            LeftActivation = 1-RightActivation;
        end
        
            
            if (State_L(end)==N) && Lpz(end,1) > L/2 % if length of last subd > length of unit-cell, no activation of N+1 possible
                RightActivation=1;
                LeftActivation = 1-RightActivation;
            elseif (State_R(end)==N) && Lpz(end,1) > L/2 % if length of last subd > length of unit-cell, no activation of N+1 possible
                RightActivation=0;
                LeftActivation = 1-RightActivation;
            end
            
            if (State_L(end)==N+1)
                RightActivation=1;
                LeftActivation = 1-RightActivation;
            elseif (State_R(end)==N+1)
                RightActivation=0;
                LeftActivation = 1-RightActivation;
            end
            
            if LeftActivation
                %[length_act,jkl,~]= fminsearch(@find_length_act_L,L/2,opts);
                [length_act,jkl,~]= fminbnd(@find_length_act_L,0,L/2,opts);
                if length_act < L/2 && jkl < 0.001 % Left activation possible
                    State_L=[State_L (State_L(end)+1)];
                    
                    g0=[];
                    gsup=gpzMax(i,State_R(end),State_L(end));
                    if (size(gpzMax,2)~=max(State_L(end),State_R(end)))
                        [g01,jkl1,~]= fminsearch(@find_g0_act_L,ginf,opts);
                        [g0,jkl,~]=fminbnd(@find_g0_act_L,ginf,gsup,opts);
                        if g01 > ginf && g01 < gsup && abs(g01-g0)>abs(g01)/20 && jkl<5
                            [g0,jkl,~]= fminsearch(@find_g0_act_L,ginf,opts);
                        end
                    else
                        jkl=10;
                    end
                    
                    if  isempty(g0) || jkl>5
                        [g0,jkl,~]= fminsearch(@find_g0_act_L,ginf,opts);
                        if g0 < ginf || jkl>5
                            divGsup=1;
                            jkl=1;
                            while jkl>0.1
                                divGsup=divGsup*2;
                                [g0,jkl,~]= fminbnd(@find_g0_act_L,ginf,ginf*divGsup,opts);
                            end
                        end
                    end
                    
                    gR=gR_end;
                    gL=gm(State_L(end-1));
                else % deactivation in the center first
                    State_L = State_L(2:end);
                    State_R = State_R(2:end);
                    
                    if nL==1
                        Psi0L=0;
                        lL0=0;
                    else
                        Psi0L=-P0pz(i,nL-1);
                        lL0=Lpz(i,nL-1);
                    end
                    
                    if nR==1
                        Psi0R=0;
                        lR0=0;
                    else
                        Psi0R=P0pz(i,nR-1);
                        lR0=Lpz(i,nR-1);
                    end
                    
                    %[LnR,jkl,~]= fminsearch(@find_P,(L/2-(lR0+lL0))/2,opts);
                    %[LnR,jkl,~]= fminbnd(@find_P,0,(L/2-(lR0+lL0)),opts);
                    
                    %[LnR,jkl,~]= fminsearch(@find_P,(L/2-(lR0+lL0))/2,opts);
                    [LnR,jkl,~]= fminsearch(@find_P,LRA,opts);
                    if LnR < 0 || jkl>0.1
                        [LnR,jkl,~]= fminbnd(@find_P,0,(L/2-(lR0+lL0)),opts);
                    end
                    
                    gL=gL_end;
                    gR=gR_end;
                    PsiR=PR_end;
                    
                end
            elseif RightActivation
                %[length_act,jkl,~]= fminsearch(@find_length_act_R,L/2,opts);
                [length_act,jkl,~]= fminbnd(@find_length_act_R,0,L/2,opts);
                if length_act < L/2 && jkl < 0.001 % Right activation possible
                    State_R=[State_R (State_R(end)+1)];
                    
                    g0=[];
                    gsup=gpzMax(i,State_R(end),State_L(end));
                    if (size(gpzMax,2)~=max(State_L(end),State_R(end)))
                        [g01,jkl1,~]= fminsearch(@find_g0_act_R,ginf,opts);
                        [g0,jkl,~]=fminbnd(@find_g0_act_R,ginf,gsup,opts);
                        if g01 > ginf && g01 < gsup && abs(g01-g0)>abs(g01)/20
                            [g0,jkl,~]= fminsearch(@find_g0_act_R,ginf,opts);
                        end
                    else
                        jkl=10;
                    end
                    
                    
                    if  isempty(g0) || jkl>5
                        [g0,jkl,~]= fminsearch(@find_g0_act_R,ginf,opts);
                        if g0 < ginf || jkl>5
                            divGsup=1;
                            jkl=1;
                            while jkl>0.1
                                divGsup=divGsup*2;
                                [g0,jkl,~]= fminbnd(@find_g0_act_R,ginf,ginf*divGsup,opts);
                            end
                        end
                    end
                    
                    gL=gL_end;
                    gR=gm(State_R(end-1));
                else % deactivation in the center first
                    State_R = State_R(2:end);
                    State_L = State_L(2:end);
                    
                    if nL==1
                        Psi0L=0;
                        lL0=0;
                    else
                        Psi0L=-P0pz(i,nL-1);
                        lL0=Lpz(i,nL-1);
                    end
                    
                    if nR==1
                        Psi0R=0;
                        lR0=0;
                    else
                        Psi0R=P0pz(i,nR-1);
                        lR0=Lpz(i,nR-1);
                    end
                    
                    %[LnR,jkl,~]= fminsearch(@find_P,(L/2-(lR0+lL0))/2,opts);
                    [LnR,jkl,~]= fminsearch(@find_P,LRA,opts);
                    if LnR < 0 || jkl>0.1
                        [LnR,jkl,~]= fminbnd(@find_P,0,(L/2-(lR0+lL0)),opts);
                    end
                    
                    gL=gL_end;
                    gR=gR_end;
                    PsiR=PR_end;
                    
                end
            end
            
        end
        
        [epsilon_inf, sigma_inf] = post_process(PsiR,gL,gR);
        SS_curve=[SS_curve;epsilon_inf,sigma_inf];
        g_crit=[g_crit g0];
        
    end
    
%     figure(8)
%     hold on
%     plot(epsilon_inf,sigma_inf,'*')
    
    
    ShearStressField();
    
    if sigma_inf<max(SS_curve(:,2))/100
        break
    end
    
end

SS_steps_R=SS_steps;
SS_steps(end,1)=SS_curve(end,1);

% figure(8)
% hold on
% plot(SS_steps(:,1),SS_steps(:,2),'-.')
% hold on
% plot(SS_steps(1:ng:end,1),SS_steps(1:ng:end,2),'o')

X=SS_curve(:,1);
Y=SS_curve(:,2);
[strength,ind]=max(Y);
eps_softening=X(ind);
if max(X)>2.1
    eps_softening=max(X);
end

% 3D PLOTS ----------------------------------------------------------------
%SurfacePlot( S3Dplot )
%SurfacePlot( P3Dplot )
%SurfacePlot( Xa3Dplot )
%SurfacePlot( Xb3Dplot )

return
%[ SS_curve,strength,eps_softening,vec_l_int ];


















    function [epsilon_inf, sigma_inf] = post_process(PsiR,gL,gR)
        
        Xinf=PsiR/(2*ET/(ET+DT)-K);
        epsilon_inf=100*(tm/(L)*((gL*(EE+DE)*(ET+DT)+gR*(EE-DE)*(ET-DT))/(EE*ET+DE*DT))+2*Xinf*ET/(EE*ET+DE*DT));
        sigma_inf=Xinf;%*ET/(ET+tm);
        
    end

    function [] = ShearStressField()

        
        %for gj=g_crit(end-1):(g_crit(end)-g_crit(end-1))/ng:g_crit(end) % evolution of g between the two last critical values
        
        %for gj=exp(log(g_crit(end-1))+(log(g_crit(end))-log(g_crit(end-1)))/ng:(log(g_crit(end))-log(g_crit(end-1)))/ng:log(g_crit(end))-(log(g_crit(end))-log(g_crit(end-1)))/ng) % evolution of g between the two last critical values
        for gj=exp(log(g_crit(end-1)):(log(g_crit(end))-log(g_crit(end-1)))/ng:log(g_crit(end))) % evolution of g between the two last critical values
            xFieldR=[];
            ShearStressR=[];
            PsiStressR=[];
            
            g0iR=gj;
            S0iR=interp1([0 gm'],[0 Sm'],g0iR,'linear');
            P0iR=0;
            sumPreviousR=0;
            for njR=State_R_previous(1):State_R_previous(end)-1 % loop on Right subdomains
                
                [liR,~]=LengthSing(type(njR),Y(njR),g0iR,gm(njR),S0iR,...
                    Sm(njR),P0iR,ET,DT,tm,EE,DE);
                
                for lengthjR=0:liR/nx:liR % loop for different values of x along the geometry
                    [gR_end,SR_end,PsiR]=field_2Stiffness_nonIterative(0,njR,EE,DE,0,0,0,Y,0,ET,DT,tm,g0iR,S0iR,P0iR,0,lengthjR,type(njR));
                    xFieldR=[xFieldR sumPreviousR+lengthjR];
                    ShearStressR=[ShearStressR SR_end];
                    PsiStressR=[PsiStressR PsiR];
                end
                if ~isempty(xFieldR)
                    sumPreviousR=xFieldR(end);
                    P0iR=PsiR;
                    g0iR=gm(njR);
                    S0iR=Sm(njR);
                end
            end
            
            xFieldL=[];
            ShearStressL=[];
            PsiStressL=[];
            
            g0iL=gj;
            S0iL=interp1([0 gm'],[0 Sm'],g0iL,'linear');
            P0iL=0;
            sumPreviousL=0;
            for njL=State_L_previous(1):State_L_previous(end)-1 % loop on Left subdomains
                
                [liL,~]=LengthSing(type(njL),Y(njL),g0iL,gm(njL),S0iL,...
                    Sm(njL),abs(P0iL),ET,DT,tm,EE,DE);
                
                for lengthjL=0:liL/nx:liL % loop for different values of x along the geometry
                    [gL_end,SL_end,PsiL]=field_2Stiffness_nonIterative(0,njL,EE,DE,0,0,0,Y,0,ET,DT,tm,g0iL,S0iL,P0iL,0,-lengthjL,type(njL));
                    xFieldL=[xFieldL sumPreviousL+lengthjL];
                    ShearStressL=[ShearStressL SL_end];
                    PsiStressL=[PsiStressL PsiL];
                end
                if ~isempty(xFieldL)
                    sumPreviousL=xFieldL(end);
                    P0iL=PsiL;
                    g0iL=gm(njL);
                    S0iL=Sm(njL);
                end
            end
            
            lL0Plot=sumPreviousL;
            lR0Plot=sumPreviousR;
            P0iLPlot=P0iL;
            P0iRPlot=P0iR;
            g0iLPlot=g0iL;
            g0iRPlot=g0iR;
            S0iLPlot=S0iL;
            S0iRPlot=S0iR;
            
            LRB=lengthR(end);
            
            %[LnR,jkl,~]= fmincon(@find_P_fixed_g0,0,[],[],[],[],0,(L/2-(lR0Plot+lL0Plot)),[],opts);
            
            [LnR,jkl,~]= fminsearch(@find_P_fixed_g0,0,opts);
            if LnR<0
                [LnR,jkl,~]= fminsearch(@find_P_fixed_g0,(L/2-(lR0Plot+lL0Plot)),opts);
                if LnR>(L/2-(lR0Plot+lL0Plot))
                    1;
                   %[LnR,jkl,~]= fminbnd(@find_P_fixed_g0,0,(L/2-(lR0Plot+lL0Plot)),opts);
                    
                    
%                             divGsup=1;
%                             jkl=1;
%                             while jkl>0.1
%                                 divGsup=divGsup*2;
%                                 [g0,jkl,~]= fminbnd(@find_P_fixed_g0,ginf,ginf*divGsup,opts);
%                             end
                            
                end
            end
            
            if (State_R_previous_previous(end)~=State_R_previous(end))
                    LRA=0;
            end
            
            [LnR,jkl,~]= fminbnd(@find_P_fixed_g0,min(LRA,LRB),max(LRA,LRB),opts);
            
            1;
            
% [LnR,jkl,~]= fminsearch(@find_P_fixed_g0,(L/2-(lR0Plot+lL0Plot)),opts);
% if LnR>(L/2-(lR0Plot+lL0Plot))
%     [LnR,jkl,~]= fminsearch(@find_P_fixed_g0,(L/2-(lR0Plot+lL0Plot))/2,opts);
%     if LnR>(L/2-(lR0Plot+lL0Plot)) || LnR<0
%         [LnR,jkl,~]= fminsearch(@find_P_fixed_g0,0,opts);
%         if LnR<0
%             [LnR,jkl,~]= fminbnd(@find_P_fixed_g0,0,(L/2-(lR0Plot+lL0Plot)),opts);
%         end
%     end
% end


            
            njR=State_R_previous(end); % last Right subdomain
            for lengthjR=0:LnR/nx:LnR % loop for different values of x along the geometry
                [gR_end,SR_end,PsiR]=field_2Stiffness_nonIterative(0,njR,EE,DE,0,0,0,Y,0,ET,DT,tm,g0iR,S0iR,P0iR,0,lengthjR,type(njR));
                xFieldR=[xFieldR sumPreviousR+lengthjR];
                ShearStressR=[ShearStressR SR_end];
                PsiStressR=[PsiStressR PsiR];
            end
            
            
            njL=State_L_previous(end); % last Left subdomain
            if isempty(xFieldR) && isempty(xFieldL)
                LnL=L/2;
            elseif isempty(xFieldR)
                LnL=L/2-(xFieldL(end));
            elseif isempty(xFieldL)
                LnL=L/2-(xFieldR(end));
            else
                LnL=L/2-(xFieldR(end)+xFieldL(end));
            end
            for lengthjL=0:LnL/nx:LnL % loop for different values of x along the geometry
                [gL_end,SL_end,PsiL]=field_2Stiffness_nonIterative(0,njL,EE,DE,0,0,0,Y,0,ET,DT,tm,g0iL,S0iL,P0iL,0,-lengthjL,type(njL));
                xFieldL=[xFieldL sumPreviousL+lengthjL];
                ShearStressL=[ShearStressL SL_end];
                PsiStressL=[PsiStressL PsiL];
            end
            
            
            %%% PLOT SHEAR STRESS FIELD %%%
%             figure(16)
%             hold on
%             plot(max(xFieldL)-xFieldL,ShearStressL)
%             plot(max(xFieldL)+xFieldR,ShearStressR)
            
            if isempty(PsiStressR)
                PR=0;
            else
                PR=PsiStressR(end);
            end
            
            S3DplotTemp=[ones(1,length(ShearStressL)+length(ShearStressR))*gj*100 ; fliplr(max(xFieldL)-xFieldL) , max(xFieldL)+xFieldR ; fliplr(ShearStressL) , ShearStressR];
            P3DplotTemp=[ones(1,length(PsiStressL)+length(PsiStressR))*gj*100 ; fliplr(max(xFieldL)-xFieldL) , max(xFieldL)+xFieldR ; fliplr(PsiStressL) , PsiStressR];
            Xa3DplotTemp=[ones(1,length(PsiStressL)+length(PsiStressR))*gj*100 ; fliplr(max(xFieldL)-xFieldL) , max(xFieldL)+xFieldR ; -([fliplr(PsiStressL) PsiStressR]+PR*((K-ET/((ET+DT)/2))/(-K+ET/((ET+DT)/2))))/(1+((ET-DT)/2)/((ET+DT)/2))];
            Xb3DplotTemp=[ones(1,length(PsiStressL)+length(PsiStressR))*gj*100 ; fliplr(max(xFieldL)-xFieldL) , max(xFieldL)+xFieldR ; ([fliplr(PsiStressL) PsiStressR]+PR*((K+ET/((ET-DT)/2))/(-K+ET/((ET+DT)/2))))/(1+((ET+DT)/2)/((ET-DT)/2))];
            
            
            S3Dplot=[S3Dplot;S3DplotTemp'];
            P3Dplot=[P3Dplot;P3DplotTemp'];
            Xa3Dplot=[Xa3Dplot;Xa3DplotTemp'];
            Xb3Dplot=[Xb3Dplot;Xb3DplotTemp'];
            
            [epsilon_step, sigma_step] = post_process(PsiR,gL_end,gR_end);
            SS_steps=[SS_steps;epsilon_step, sigma_step];
            
        end
        
    end


    function find_P_fixed_g0 = find_P_fixed_g0(LnR_guess_fixed_g0)

        LnL_guess_fixed_g0=L/2-(lR0Plot+lL0Plot+LnR_guess_fixed_g0);
        
        [~,~,PL_end_plot]=field_2Stiffness_nonIterative(State_L_previous,State_L_previous(end),EE,DE,G,Sm,gm,Y,0,ET,DT,tm,...
            g0iLPlot,S0iLPlot,P0iLPlot,0,-LnL_guess_fixed_g0,type(State_L_previous(end)));
        
        [~,~,PR_end_plot]=field_2Stiffness_nonIterative(State_R_previous,State_R_previous(end),EE,DE,G,Sm,gm,Y,0,ET,DT,tm,...
            g0iRPlot,S0iRPlot,P0iRPlot,0,LnR_guess_fixed_g0,type(State_R_previous(end)));
        
        %PL_end_plot
        find_P_fixed_g0=abs(PL_end_plot+p*PR_end_plot);
    end




    function error = find_g0_1sub_end(g0_guess)
        
        % Left %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lL=0;
        lengthL=[];
        P0iL=0;
        g0iL=g0_guess;
        S0iL=interp1([0 gm'],[0 Sm'],g0iL,'linear');
        for iL=0:nL
            [liL,P0iL]=LengthSing(type(i+iL),Y(i+iL),g0iL,gm(i+iL),S0iL,...
                Sm(i+iL),abs(P0iL),ET,DT,tm,EE,DE);
            g0iL=gm(i+iL);
            S0iL=Sm(i+iL);
            lL=lL+liL;
            lengthL=[lengthL liL];
        end
        PsiL=-P0iL;
        
        % Right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lR=0;
        lengthR=[];
        P0iR=0;
        g0iR=g0_guess;
        S0iR=interp1([0 gm'],[0 Sm'],g0iR,'linear');
        for iR=0:nR
            [liR,P0iR]=LengthSing(type(i+iR),Y(i+iR),g0iR,gm(i+iR),S0iR,...
                Sm(i+iR),P0iR,ET,DT,tm,EE,DE);
            g0iR=gm(i+iR);
            S0iR=Sm(i+iR);
            lR=lR+liR;
            lengthR=[lengthR liR];
        end
        PsiR=P0iR;
        
        error=abs(lR+lL-L/2);
    end



    function error = find_g0_act_L(g0_guess)
        
        % Left %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lL=0;
        lengthL=[];
        P0iL=0;
        g0iL=g0_guess;
        S0iL=interp1([0 gm'],[0 Sm'],g0iL,'linear');
        for iL=0:nL
            [liL,P0iL]=LengthSing(type(i+iL),Y(i+iL),g0iL,gm(i+iL),S0iL,...
                Sm(i+iL),abs(P0iL),ET,DT,tm,EE,DE);
            g0iL=gm(i+iL);
            S0iL=Sm(i+iL);
            lL=lL+liL;
            lengthL=[lengthL liL];
        end
        PsiL=-P0iL;
        
        % Right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lR=0;
        lengthR=[];
        P0iR=0;
        g0iR=g0_guess;
        S0iR=interp1([0 gm'],[0 Sm'],g0iR,'linear');
        for iR=0:nR-1
            [liR,P0iR]=LengthSing(type(i+iR),Y(i+iR),g0iR,gm(i+iR),S0iR,...
                Sm(i+iR),P0iR,ET,DT,tm,EE,DE);
            g0iR=gm(i+iR);
            S0iR=Sm(i+iR);
            lR=lR+liR;
            lengthR=[lengthR liR];
        end
        [gR_end,SR_end,PsiR]=field_2Stiffness_nonIterative(State_R,i+nR,EE,DE,G,Sm,gm,Y,0,ET,DT,tm,g0iR,S0iR,P0iR,0,L/2-(lL+lR),type(i+nR));
        lengthR=[lengthR L/2-(lL+lR)];
        
        error=abs(PsiL+p*PsiR)+g0_guess/1000;
    end

    function error = find_g0_act_R(g0_guess)
        
        % Right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lR=0;
        lengthR=[];
        P0iR=0;
        g0iR=g0_guess;
        S0iR=interp1([0 gm'],[0 Sm'],g0iR,'linear');
        for iR=0:nR
            [liR,P0iR]=LengthSing(type(i+iR),Y(i+iR),g0iR,gm(i+iR),S0iR,...
                Sm(i+iR),P0iR,ET,DT,tm,EE,DE);
            g0iR=gm(i+iR);
            S0iR=Sm(i+iR);
            lR=lR+liR;
            lengthR=[lengthR liR];
        end
        PsiR=P0iR;
        
        % Left %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lL=0;
        lengthL=[];
        P0iL=0;
        g0iL=g0_guess;
        S0iL=interp1([0 gm'],[0 Sm'],g0iL,'linear');
        for iL=0:nL-1
            [liL,P0iL]=LengthSing(type(i+iL),Y(i+iL),g0iL,gm(i+iL),S0iL,...
                Sm(i+iL),abs(P0iL),ET,DT,tm,EE,DE);
            g0iL=gm(i+iL);
            S0iL=Sm(i+iL);
            lL=lL+liL;
            lengthL=[lengthL liL];
        end
        [gL_end,SL_end,PsiL]=field_2Stiffness_nonIterative(State_L,i+nL,EE,DE,G,Sm,gm,Y,0,ET,DT,tm,g0iL,S0iL,-P0iL,0,-(L/2-(lL+lR)),type(i+nL));
        lengthL=[lengthL L/2-(lL+lR)];
        
        error=abs(PsiL+p*PsiR)+g0_guess/1000;
    end



    function error = find_length_act_L(length_guess)
        
        % Left %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lL=0;
        lengthL=[];
        P0iL=0;
        for iL=0:nL-1
            [liL,P0iL]=LengthSing(type(i+iL+1),Y(i+iL+1),gm(i+iL),gm(i+iL+1),Sm(i+iL),...
                Sm(i+iL+1),abs(P0iL),ET,DT,tm,EE,DE);
            lL=lL+liL;
            lengthL=[lengthL liL];
        end
        PsiL=-P0iL;
        
        % Right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lR=0;
        lengthR=[];
        P0iR=0;
        for iR=0:nR-1-1
            [liR,P0iR]=LengthSing(type(i+iR+1),Y(i+iR+1),gm(i+iR),gm(i+iR+1),Sm(i+iR),...
                Sm(i+iR+1),P0iR,ET,DT,tm,EE,DE);
            lR=lR+liR;
            lengthR=[lengthR liR];
        end
        [gR_end,SR_end,PsiR]=field_2Stiffness_nonIterative(State_R,i+nR,EE,DE,G,Sm,gm,Y,0,ET,DT,tm,gm(i+nR-1),Sm(i+nR-1),P0iR,0,length_guess-(lL+lR),type(i+nR));
        lengthR=[lengthR length_guess-(lL+lR)];
        
        error=abs(PsiL+p*PsiR);
    end



    function error = find_length_act_R(length_guess)
        
        % Left %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lR=0;
        lengthR=[];
        P0iR=0;
        for iR=0:nR-1
            [liR,P0iR]=LengthSing(type(i+iR+1),Y(i+iR+1),gm(i+iR),gm(i+iR+1),Sm(i+iR),...
                Sm(i+iR+1),P0iR,ET,DT,tm,EE,DE);
            lR=lR+liR;
            lengthR=[lengthR liR];
        end
        PsiR=P0iR;
        
        % Right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lL=0;
        lengthL=[];
        P0iL=0;
        for iL=0:nL-1-1
            [liL,P0iL]=LengthSing(type(i+iL+1),Y(i+iL+1),gm(i+iL),gm(i+iL+1),Sm(i+iL),...
                Sm(i+iL+1),abs(P0iL),ET,DT,tm,EE,DE);
            lL=lL+liL;
            lengthL=[lengthL liL];
        end
        [gL_end,SL_end,PsiL]=field_2Stiffness_nonIterative(State_L,i+nL,EE,DE,G,Sm,gm,Y,0,ET,DT,tm,gm(i+nL-1),Sm(i+nL-1),-P0iL,0,-(length_guess-(lL+lR)),type(i+nL));
        lengthL=[lengthL length_guess-(lL+lR)];
        
        error=abs(PsiL+p*PsiR);
    end

    function [] = find_Psi_end(g0_known)
        
        % Right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lR=0;
        lengthR=[];
        P0iR=0;
        g0iR=g0_known;
        S0iR=interp1([0 gm'],[0 Sm'],g0iR,'linear');
        for iR=0:nR
            [liR,P0iR]=LengthSing(type(i+iR),Y(i+iR),g0iR,gm(i+iR),S0iR,...
                Sm(i+iR),P0iR,ET,DT,tm,EE,DE);
            g0iR=gm(i+iR);
            S0iR=Sm(i+iR);
            lR=lR+liR;
            lengthR=[lengthR liR];
        end
        PsiR=P0iR;
        
        % Left %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lL=0;
        lengthL=[];
        P0iL=0;
        g0iL=g0_known;
        S0iL=interp1([0 gm'],[0 Sm'],g0iL,'linear');
        for iL=0:nL
            [liL,P0iL]=LengthSing(type(i+iL),Y(i+iL),g0iL,gm(i+iL),S0iL,...
                Sm(i+iL),abs(P0iL),ET,DT,tm,EE,DE);
            g0iL=gm(i+iL);
            S0iL=Sm(i+iL);
            lL=lL+liL;
            lengthL=[lengthL liL];
        end
        PsiL=-P0iL;
        
    end






    function find_P = find_P(LnR_guess)

        LnL_guess=L/2-(lR0+lL0+LnR_guess);
        [gL_end,SL_end,PL_end]=field_2Stiffness_nonIterative(State_L,State_L(end),EE,DE,G,Sm,gm,Y,0,ET,DT,tm,gm(State_L(end)-1),Sm(State_L(end)-1),Psi0L,0,-LnL_guess,type(State_L(end)));
        [gR_end,SR_end,PR_end]=field_2Stiffness_nonIterative(State_R,State_R(end),EE,DE,G,Sm,gm,Y,0,ET,DT,tm,gm(State_R(end)-1),Sm(State_R(end)-1),Psi0R,0,LnR_guess,type(State_R(end)));
        
        lengthR=LnR_guess;
        lengthL=LnL_guess;
        
        find_P=abs(PL_end+p*PR_end);
    end


%     function find_DXL = find_DXL(guess)
%         Xinf_guess=guess(1);
%         l_int_guess=max(0,guess(2));
%         l_int_guess=min(L/2,l_int_guess);
%
%         DX0_int_R=Xinf_guess*2*(DE*ET)/(EE*ET+DE*DT);
%         DX0_int_L=-DX0_int_R;
%
%         %f_int=@(l_int_guess)find_l_int(l_int_guess,Xinf_guess,DX0_int_R,DX0_int_L);
%         %[l_int,fval,exitflag]=fminbnd(f_int,l_char_ini-0.1,l_char_ini+0.1,opts); % Use of function find_l_int
%
%         [DX_L,X_L,g_L,S_L,Dl_L,l_L,Dli_L,sol_is_zero_L]=resolution_2Stiffness(l_int_guess,State_L,gini,Xinf_guess,DX0_int_L,opts,ET,-DT,tm,EE,-DE,G,Sm,gm,Y);
%         [DX_R,X_R,g_R,S_R,Dl_R,l_R,Dli_R,sol_is_zero_R]=resolution_2Stiffness(L/2-l_int_guess,State_R,gini,Xinf_guess,DX0_int_R,opts,ET,DT,tm,EE,DE,G,Sm,gm,Y);
%
%         find_DXL=abs(DX_R(end)*(ET+DT)/2-Xinf_guess*ET)+abs(DX_L(end)*(ET-DT)/2-Xinf_guess*ET);
%     end

%     function find_l_int = find_l_int(l_int_guess,Xinf_guess,DX0_int_R,DX0_int_L)
%         [DX_L,X_L,g_L,S_L,Dl_L,l_L,Dli_L,sol_is_zero_L]=resolution_2Stiffness(l_int_guess,State_L,gini,Xinf_guess,DX0_int_L,opts,ET,-DT,tm,EE,-DE,G,Sm,gm,Y);
%         [DX_R,X_R,g_R,S_R,Dl_R,l_R,Dli_R,sol_is_zero_R]=resolution_2Stiffness(L/2-l_int_guess,State_R,gini,Xinf_guess,DX0_int_R,opts,ET,DT,tm,EE,DE,G,Sm,gm,Y);
%
%         find_l_int=abs((ET-DT)*DX_L(end)-(ET+DT)*DX_R(end));
%         %find_l_int=(DX_L(end)-DX_R(end));
%     end

% Add last point to curve, where Xinf = 0
if SS_curve(end,1) <0
    SS_curve=SS_curve(1:end-1,:);
end
if abs(SS_curve(end,2)) > max(SS_curve(:,2))/10
    %%% ENERGY
    Nfin=10;
    [~,~,~,~,~,~,~,~,~,Gic] = inputs_2Stiffness;
    de=SS_curve(2:end,1)-SS_curve(1:end-1,1);
    se=(SS_curve(2:end,2)+SS_curve(1:end-1,2))/2;
    GicTemp=sum(de.*se);
    dG=Gic/(ET+tm)*100-GicTemp;
    depsilon=[2*dG/SS_curve(end,2)*[1/Nfin:1/Nfin:1]]';
    dsigma=[SS_curve(end,2)*[1/Nfin:1/Nfin:1]]';
    SS_curve=[SS_curve ; SS_curve(end,1)+depsilon SS_curve(end,2)-dsigma];
else
    %%% LINEAR
    SS_curve=SS_curve(1:end-1,:);
    XF=SS_curve(end,1)-(SS_curve(end,1)-SS_curve(end-1,1))*(SS_curve(end,2)/(SS_curve(end,2)-SS_curve(end-1,2)));
    SS_curve=[SS_curve ; XF 0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%  ONSET OF SOFTENING:  %
%%%%%%%%%%%%%%%%%%%%%%%%%

X=SS_curve(:,1);
Y=SS_curve(:,2);
[strength,ind]=max(Y);
eps_softening=X(ind);
if max(X)>2.1
    eps_softening=max(X);
end

% PLOT ----------------
% figure(8);
% hold on
% plot(SS_curve(:,1),SS_curve(:,2),'go')


end




%%%%%%%%%%%%%%%%%%%%%%%%%
%  ONSET OF SOFTENING:  %
%%%%%%%%%%%%%%%%%%%%%%%%%

% [SS,lint]=general_loop_2Stiffness;
% X=SS(:,1);
% Y=SS(:,2);
% [S,ind]=max(Y);
% e=X(ind);
% [e S]

%%%%%%%%%%%%%%%%%%%%%%%%%
%  plot versus DT--DE:  %
%%%%%%%%%%%%%%%%%%%%%%%%%

% I=[];
% S=[];
% e=[];
% lim=.9;
% N=10;
% AE=0;
% for i=-lim:lim/N:lim
% SS=general_loop_2Stiffness_nonIterative(AE,i);
% XSS=SS(:,1);
% YSS=SS(:,2);
% [S_t,ind]=max(YSS);
% e_t=XSS(ind);
% if max(XSS)>XSS(end)
% e_t=max(XSS);
% end
% I=[I i];
% S=[S S_t];
% e=[e e_t];
% end
% figure
% plot(I,e)
% figure
% plot(I,S)