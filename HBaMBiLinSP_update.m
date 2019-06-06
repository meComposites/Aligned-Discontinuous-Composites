function [stateEnd,deInt,Epsilon,Sigma] = ...
    HBaMBiLinSP_update(SumE,DeltaE,SumT,DeltaT,tI,Lo,SI,GiicI,GI,DeltaSigma_min)
% dbstop in file
%FE in paper (Fig.5a)
%HBaMBiLinSP(200000,40000,0.2,0.06,0.01,2.0,80,2.0,1500,100);

%non hybrid, my paper:
%HBaMBiLinSP(200000,0.01,0.2,0,0.01,2.0,50,1,1000,100)

%non hybrid, my paper:
%HBaMBiLinSP(200000,0.01,0.2,0,0.01,10.0,50,1,1000,100)


%%
DeltaT=-DeltaT;
DeltaE=-DeltaE;

TA=(SumT-DeltaT)/2;
TB=(SumT+DeltaT)/2;
EA=(SumE-DeltaE)/2;
EB=(SumE+DeltaE)/2;

KA=TA*EA;
KB=TB*EB;

eta=KA/KB;  % Eq.1
Teq=((SumT.^2-DeltaT.^2)./(2*SumT)); % Equivalent thickness

gI1=SI/GI;
gI2=2*GiicI/(tI*SI);

Keq=2*(KA*KB)/(KA+KB);
YI1=sqrt(2*GI/(tI*Keq)); %Eq.2
YI2=sqrt(2*SI/(gI2-gI1)/(tI*Keq)); %Eq.2

kappa=2*(DeltaE*SumT)/(KA+KB); %Eq.2


%% Linear-elastic domain: state=[1]
%  (Based on Stiffness paper)

LR_1=1/(2*YI1)*log((1/eta+exp(YI1*Lo))/(1/eta+exp(-YI1*Lo))); %Eq.22 
    %(replacing LL->LR, eta->1/eta, 2*L->Lo)

% Defining the transition point from [1] to [1 2]: tauR=SI or gR=gI1 
tauR(1)=SI;
tau0_1_12=SI/cosh(YI1*LR_1);  %Eq.9b-c

% Defining tauL, PsiR for the same tau0_1_12
PsiR(1)=2/(YI1*Teq)*sqrt(tauR(1)^2-tau0_1_12^2); % SS-paper, Tab.2 Eq.a
tauL(1)=tauR(1)/eta*sqrt(1+(eta^2-1)/cosh(YI1*LR_1)^2); % calculated
    %from BC Eq.6, defining PsiL and PsiR from Tab.2 Eq.a, tau0 as above

% Defining required accuracy for numerical solutions
opts=optimset('TolX',10^(-18)*min(tau0_1_12,1),'TolFun',10^(-18)*Lo);

%% Calculating evolution of status. Starts by [1]->[1 2]->?
% Full PZ length for [2] on the right: tau0=SI, Psi0=0, tauR=0 in Tab.3, Eq.c
L_2R_full=pi/(YI2*2);
L_2L_when2Rfull=2/YI2*atan(sqrt((1-sqrt(1-1/eta^2))/(1+sqrt(1-1/eta^2)))); 
        %Found tauL/SI=sqrt(1-1/eta^2) as follows: 
        %PsiR from Tab.2, Eq.c; PsiL from Eq.5b;
        %PsiL(LL) from Tab.1, Eq.f; LL from PsiL=PsiL(LL);
        %tauL(LL) from Tab.1, Eq.f; simplified cos(asin(z))=sqrt(1-z^2)
L_2full=L_2R_full+L_2L_when2Rfull;

%% Calculating variables to define whether [2 1 2]->[2 1 2 3] or [2]
% Assuming evolution of status is [1]->[1 2]->[2 1 2]
% Length of PZ for fully developed [2R], and equilibrium [2L] (i.e. g0>gI1)
% Tab.3, Eq.c, simplified for tau0=SI, Psi0=0)

%%
%criating guessed variables (to share with nested functions)
L_1_guess=NaN;
PsiR_guess=NaN;

%% Defining the transition point from state=[1 2] to [2 1 2] or [1 2 3]: 
% Assuming transition to [2 1 2] (with tauL=SI), and checking whether 
% tauR>0 (confirming [2 1 2]) or tauR<0 (thus transition to [1 2 3])
% (this will be stored at the end of the [1 2] interval, at 1+ni points

% initiating variables for outputs of @transition_12_212:
Psi_1_2_guess=NaN;
PsiR_guess=NaN;
tauR_guess=NaN;
L_2R_guess=NaN;

try
% calculating tau0 for transtion between [1 2] and [2 1 2]:
tau0_End12=abs(fzero(@transition_12_212,tau0_1_12*SI/tauL(1),opts));
%tau0_End12=fzero(@transition_12_212,[tau0_1_12,SI],opts);
% also calculated PsiR=PsiR_guess, tauR_guess, and L_2R=L_2R_guess
catch
    1; 
end

if tauR_guess>0
    %If tauR>0, then [1 2]->[2 1 2] (and tauL_End12=SI). In this case,
    %the variables calculated using @transition_12_212 are correct.
    stateEnd=212;
    tauR_End12=tauR_guess;
    tauL_End12=SI; % defines last point

else
    %If tauR=0 (Im solution for tau0_End12), then [1 2]->[1 2 3] (and tauR_End12=0). 
    %Correcting tauR_End12:
    tauR_End12=0;
    
    %In this case, must re-calculate transition 
    % initiating variables for outputs of @transition_12_123:
    L_1L_guess=NaN;
    PsiL_guess=NaN;
    
    try
    % calculating tau0 for transtion between [1 2] and [1 2 3]:
    tau0_End12=fzero(@transition_12_123,[tau0_1_12,SI],opts);
    catch
        1;
    end
    
    % also calculated L_1L=L_1L_guess. Thus calculating
    tauL_End12=tau0_End12*cosh(YI1*L_1L_guess);
    
    stateEnd=123;
end 

%From transition function: also calculated PsiR_End12
PsiR_End12=PsiR_guess;

%% Defining number of points to calculate for [1 2]

DeltaPsiR_12 = PsiR_End12-PsiR(1);
DeltaSigmaO_12 = DeltaPsiR_12/(SumT/TB-kappa); %Eq.7
n_12 = ceil(DeltaSigmaO_12/DeltaSigma_min);

% Assigning variables to last point for s=[1 2] 
PsiR(1+n_12,1)=PsiR_End12;
tauL(1+n_12,1)=tauL_End12;
tauR(1+n_12,1)=tauR_End12;

%% Calculating points for state=[1 2]
% calculate ni intermediate SS point:
tau0_int=logspace(log10(tau0_1_12),log10(tau0_End12),n_12+1)';
tau0_int=tau0_int(2:end-1);

tauL_int_guess=linspace(tauL(1),tauL_End12,n_12+1);
tauL_int_guess=tauL_int_guess(2:end-1);

% Calculating PsiR, gL and gR for s=[1 2]

for k=1:n_12-1
    
    tau0_now=tau0_int(k);

    % Calculating [1] on the right side
    Psi_1_2_now=Psi_1_from_tauL(YI1,Teq,tau0_now,SI);
    L_1R_now=L_1_from_tau(YI1,tau0_now,SI);
    
    try
    % Calculating PsiR and gL
%    tauL(k+1)=fzero(@within_12,[max(tauL(k),tau0_now),SI],opts);
    tauL(k+1)=fzero(@within_12,tauL_int_guess(k),opts);
    catch 
        1;
    end
    % Max required because, for very short overlaps, tauL->tau0,
    % which may lead to tau0(k+1)>tauL(k)
    PsiR(k+1)=PsiR_guess;   
    
    tauR(k+1)=tau_2_from_x(YI2,Teq,SI,Psi_1_2_now,L_2R_guess);
end

gL=tauL/GI;
gR=g_2_from_tau(gI1,gI2,SI,tauR);    

%% If stateEnd=[2 1 2], defining end of state=[2 1 2]
if stateEnd==212

    % Calculating whether we will have [2 1 2]->[2] or [2 1 2]->[2 1 2 3]
    % this depends on whether there the overlap is long enough to 
    % develop subdomain [2] completely

    % tauL          <-[2]--> tau0=SI <-[2]->     tauR=0   
    % x=-LL         <-L_2L-> x=0     <-L_2full-> x=LR
    % Psi=-PsiR/eta          Psi=0               Psi=PsiR(Psi=0,tau0=SI,tauR=0)

    if L_2full > Lo
        %% [2 1 2]->[2] (too short for subdomain [2] to fully develop)
        stateEnd=2;
        tau0_End=SI;
        
        % Defining transition point from [2 1 2] to [2]: tau0=SI, Psi_1_2=0
        Psi_1_2_now=0;
        L_1_now=0;
        
        %Lower-bound for tauL: Tab.1, Eq.f with Psi0=0, tau0=SI, x=L_2fullL
        % (i.e. when [2R] is fully developed, and 0<tauL<SI)
        tauL_min=SI*cos(YI2*L_2L_when2Rfull);
        
        try 
        % Calculating tauL_End212
        tauL_End212=fzero(@within_212,[tauL_min,SI],opts);
        catch
            1;
        end
        % also calculated PsiR=PsiR_guess and L_2R=L_2R_guess
        PsiR_End212=PsiR_guess;
        tauR_End212=tau_2_from_x(YI2,Teq,SI,Psi_1_2_now,L_2R_guess);
        
    else
        %% [2 1 2]->[2 1 2 3] (subdomain [2] can develop fully)
        stateEnd=2123;
        
        %Defining transition point from [2 1 2] to [2 1 2 3]: tauR=0
        tauR_End212=0; % defines last point
        
        % find tau0_212_2123 (which will be the maximum value for tau0)
        L_2L_guess=0;
        
        try
        tau0_End=fzero(@transition_212_2123,[tau0_End12,SI],opts); 
        catch 
            1;
        end
        % also calculated PsiR=PsiR_guess and L_2R=L_2R_guess
        PsiR_End212=PsiR_guess;
        tauL_End212=tau_2_from_x(YI2,Teq,SI,Psi_1_2_guess,L_2L_guess);        

    end % end of checking how [2 1 2] will end
    
    %% Defining number of points to calculate for [2 1 2]

    DeltaPsiR_212 = PsiR_End212-PsiR(1);
    DeltaSigmaO_212 = DeltaPsiR_212/(SumT/TB-kappa); %Eq.7
    n_212 = ceil(DeltaSigmaO_212/DeltaSigma_min);

    % Assigning variables to last point for s=[2 1 2] 
    PsiR(1+n_12+n_212)=PsiR_End212;
    tauL(1+n_12+n_212)=tauL_End212;
    tauR(1+n_12+n_212)=tauR_End212;
    
    %% Calculating points for s=[2 1 2] (for both options Lo_2full<>Lo)

    % tau0 between tau0_12_212 (onset [2 1 2]) and tau0_End (end [2 1 2])
    tau0_int=logspace(log10(tau0_End12),log10(tau0_End),n_212+1);
    tau0_int=tau0_int(2:end-1);
    
    tauL_int_guess=linspace(tauL_End212,SI,n_212+1);
    tauL_int_guess=tauL_int_guess(2:end-1);

    
    for k=1:length(tau0_int)
        tau0_now=tau0_int(k);
        % Calculating central subdomain (from tau=tau0 to tau=SI)
        Psi_1_2_now=Psi_1_from_tauL(YI1,Teq,tau0_now,SI);
        L_1_now=L_1_from_tau(YI1,tau0_now,SI);
        
        % Calculating tauL, PsiR, tauR
        L_2R_guess=0;
        
        try
        %tauL(n_12+k+1)=fzero(@within_212,[tauL_End212,SI],opts);
        tauL(n_12+k+1)=fzero(@within_212,tauL_int_guess(k),opts);
        catch 
            1;
        end
        PsiR(n_12+k+1)=PsiR_guess;
        
        tauR(n_12+k+1)=tau_2_from_x(YI2,Teq,SI,Psi_1_2_now,L_2R_guess);
    end
    
    %Calculating the shear strains for [2 1 2]
    gL(2+n_12:1+n_12+n_212)=g_2_from_tau(gI1,gI2,SI,tauL(2+n_12:end));
    gR(2+n_12:1+n_12+n_212)=g_2_from_tau(gI1,gI2,SI,tauR(2+n_12:end));
end

%% Post-processing
[Epsilon,Sigma]=post_process(PsiR,gL,gR,kappa,KA,KB,TB,SumT,tI,Lo);
deInt=gI2*tI/Lo-Epsilon(end);
%%%% SORAIA - modified to calculate deInt analytically. Small
%%%% difference, but this one is the most correct one!

return;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identification of functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
    function BCerror_12_212 = transition_12_212(tau0_guess)
        % tauL=SI       <-[1]-> tau0  <-[1]-> SI           <-[2]--> tauR   
        % x=-L_1        <-L_1-> x=0   <-L_1-> x=L_1        <-L_2R-> x=LR
        % Psi=-Psi_1_2          Psi=0         Psi=Psi_1_2           Psi=PsiR
        
        % Calculating [1] on the left side (from tau=tau0 to tau=SI):
        Psi_1_2_guess=Psi_1_from_tauL(YI1,Teq,tau0_guess,SI);
        L_1_guess=real(L_1_from_tau(YI1,tau0_guess,SI));
        
        % Calculating [1 2] on the right side
        PsiR_guess=Psi_1_2_guess*eta; %Eq.6
        tauR_guess=real(tau_2_from_Psi(YI2,Teq,SI,Psi_1_2_guess,PsiR_guess));
        L_2R_guess=real(L_2_from_tau(YI2,Teq,SI,tauR_guess,Psi_1_2_guess));
        % Re(tauR) because tauR may be Im, depending on tau0_guess 
        
        % Impose correct length for overlap
        BCerror_12_212=real(Lo-2*L_1_guess-L_2R_guess);
    end

%%
    function BCerror_12_123 = transition_12_123(tau0_guess)
        % tauL=SI   <-[1]-> tau0  <-[1]-> SI          <-[2]->      tauR=0   
        % x=-L_1L   <-L_1-> x=0   <-L_1-> x=L_1R      <-L_2fullR-> x=LR
        % Psi=PsiL          Psi=0         Psi=Psi_1_2              Psi=PsiR
        
        % Calculating [1R] (from tau=tau0 to tau=SI):
        Psi_1_2_guess=Psi_1_from_tauL(YI1,Teq,tau0_guess,SI);
        L_1R_guess=L_1_from_tau(YI1,tau0_guess,SI);
        
        % Calculating [2R] (from tau=SI to tauR=0):
        L_2R_guess=L_2_from_tau(YI2,Teq,SI,0,Psi_1_2_guess);
                    % length of subdomain [2] on the right
        PsiR_guess=Psi_2_from_L(YI2,Teq,SI,Psi_1_2_guess,L_2R_guess);

        % Calculating [1L]
        PsiL_guess=-PsiR_guess/eta; %Eq.6     
        tauL_guess=tauL_1_from_PsiL(YI1,Teq,tau0_guess,-PsiL_guess);
        L_1L_guess=L_1_from_tau(YI1,tau0_guess,tauL_guess);
        
        % Impose correct length for overlap
        BCerror_12_123=real(Lo-L_1L_guess-L_1R_guess-L_2R_guess);
    end

%%
    function BCerror_212_2123 = transition_212_2123(tau0_guess)
        % tauL=SI   <-[2]--> SI            <-[1]-> tau0  <-[1]-> SI           <-[2]--> tauR   
        % x=-LL     <-L_2L-> x=-L_1        <-L_1-> x=0   <-L_1-> x=L_1        <-L_2R-> x=LR
        % Psi=PsiL           Psi=-Psi_1_2          Psi=0         Psi=Psi_1_2           Psi=PsiR
        
        % Calculating [1] on both sides (from tau=tau0 to tau=SI):
        Psi_1_2_guess=Psi_1_from_tauL(YI1,Teq,tau0_guess,SI);
        L_1_guess=L_1_from_tau(YI1,tau0_guess,SI);
        
        % Calculating [2] on the right side (from tau=SI to tau=0):
        L_2R_guess=L_2_from_tau(YI2,Teq,SI,0,Psi_1_2_guess); 
            % length of subdomain [2] on the right 
        PsiR_guess=Psi_2_from_L(YI2,Teq,SI,Psi_1_2_guess,L_2R_guess);

        % Calculating [2] on the left side
        PsiL_guess=-PsiR_guess/eta; %Eq.6     
        tauL_guess=tau_2_from_Psi(YI2,Teq,SI,Psi_1_2_guess,-PsiL_guess);
        L_2L_guess=L_2_from_tau(YI2,Teq,SI,tauL_guess,Psi_1_2_guess);

        % Impose correct length for overlap       
        BCerror_212_2123=real(Lo-L_2L_guess-2*L_1_guess-L_2R_guess);
    end

%%
    function BCerror_12 = within_12(tauL_guess)
        % tauL     <-[1L]->  tau0  <-[1R]-> SI           <-[2]--> tauR   
        % x=-LL    <-L_1L->  x=0   <-L_1R-> x=L_1        <-L_2R-> x=LR
        % Psi=PsiL           Psi=0          Psi=Psi_1_2           Psi=PsiR
        
        % Calculating [1L] (from tau=tau0 to tau=tauL):
        PsiL_guess=-Psi_1_from_tauL(YI1,Teq,tau0_now,tauL_guess);
        L_1L_guess=L_1_from_tau(YI1,tau0_now,tauL_guess);
 
        % [1R] (from tau=tau0 to tau=SI) has been calculated already
        % Psi_1_2_now, L_1R_now
        
        % Calculating [2R]
        PsiR_guess=-PsiL_guess*eta; %Eq.6     
        tauR_guess=tau_2_from_Psi(YI2,Teq,SI,Psi_1_2_now,PsiR_guess);
        L_2R_guess=L_2_from_tau(YI2,Teq,SI,min(real(tauR_guess),SI),...
            Psi_1_2_now);
        %Min needed because, for the transition between End=2->123, 
        %this might not be verified 
        
        % Impose correct length for overlap        
        BCerror_12=real(Lo-L_1L_guess-L_1R_now-L_2R_guess);
    end
%%
    function BCerror_212 = within_212(tauL_guess)
        % tauL     <-[2]--> SI            <-[1]-> tau0  <-[1]-> SI           <-[2]--> tauR   
        % x=-LL    <-L_2L-> x=-L_1        <-L_1-> x=0   <-L_1-> x=L_1        <-L_2R-> x=LR
        % Psi=PsiL          Psi=-Psi_1_2         Psi=0         Psi=Psi_1_2          Psi=PsiR
        
        % Calculating [2] on left side (from tau=SI to tau=tauL):
        L_2L_guess=L_2_from_tau(YI2,Teq,SI,tauL_guess,Psi_1_2_now);
        PsiL_guess=-Psi_2_from_L(YI2,Teq,SI,Psi_1_2_now,L_2L_guess);

        % Calculating [2] on the right side (from x=L_1 to xLR)
        PsiR_guess=-PsiL_guess*eta; %Eq.6 
        tauR_guess=tau_2_from_Psi(YI2,Teq,SI,Psi_1_2_now,PsiR_guess);
        L_2R_guess=L_2_from_tau(YI2,Teq,SI,real(tauR_guess),Psi_1_2_now);       
            % Re(tauR) because, if tauL_guess=tauL_End212, 
            % PsiR_guess(tau0_now) is so high that tauR_guess is Im.
            % Without Re(tauR), it would only work for tau0=tau0_End212.
            
        % Impose correct length for overlap        
        BCerror_212=real(Lo-L_2L_guess-2*L_1_now-L_2R_guess);
    end


end
% end of main function
% Defining auxiliary functions (for fields) below

%% functions for fields in subdomain 1 (from tau=tau0 to tau=tauL):
%Based on SS-paper

function Psi_1=Psi_1_from_tauL(YI1,Teq,tau0,tauL)
Psi_1=2/(YI1*Teq)*sqrt(tauL^2-tau0^2); %Tab.2, Eq.a
end

function tau_1=tauL_1_from_PsiL(YI1,Teq,tau0,PsiL)
tau_1=sqrt(tau0^2+(YI1*Teq/2)^2*PsiL^2);	%Tab.2, Eq.a, 
                                            %with Psi0=0, solving wrt tauL
end

function L_1=L_1_from_tau(YI1,tau0,tauL)
L_1=2/YI1*atanh(sqrt((tauL-tau0)/(tauL+tau0))); 
    %Tab.3 Eq.a, simplified for Psi0=0
if isinf(L_1)
    L_1=1/YI1*log(2*tauL/tau0);
%L_1(L_1==Inf)=1/YI1*log(2*tauL/tau0);
end
end

%% functions for fields in subdomain 2
%Based on SS-paper

function L_2=L_2_from_tau(YI2,Teq,tau0,tauL,Psi0)
L_2=2/YI2*atan(YI2*Teq/(2*(tau0+tauL))*...
    (sqrt(Psi0^2+(2/(YI2*Teq))^2*(tau0^2-tauL^2))-Psi0)); %Tab.3 Eq.c
end

function tau_2=tau_2_from_Psi(YI2,Teq,tau0,Psi0,PsiL)
tau_2=sqrt(tau0^2-(PsiL^2-Psi0^2)*(YI2*Teq/2)^2);   %Tab.2 Eq.c,
                                                    %solved wrt tauL
end

function Psi_2=Psi_2_from_L(YI2,Teq,tau0,Psi0,x)
% all in local coordinates
Psi_2=Psi0*cos(YI2*x)+2/(YI2*Teq)*tau0*sin(YI2*x); %Tab.1, Eq.f
end

function tau_2=tau_2_from_x(YI2,Teq,tau0,Psi0,x)
tau_2=tau0*cos(YI2*x)-YI2*Teq/2*Psi0*sin(YI2*x); %Tab.1, Eq.f
end

function g_2=g_2_from_tau(gI1,gI2,SI,tau)
g_2=gI2-tau/SI*(gI2-gI1); %Tab.1, Eq.f
end

%% post-processing functions

function [epsilonO,sigmaO] = post_process(PsiR,gL,gR,...
    kappa,KA,KB,TB,SumT,tI,Lo)
sigmaO=PsiR./(SumT/TB-kappa); %Eq.7
epsilonO=tI/Lo*(gL*KB+gR*KA)./(KA+KB)+sigmaO*SumT/(KA+KB); %Eq.8
end
