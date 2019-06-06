function [AA,AAE,AAde, AB,ABE,ABde, BB,BBE,BBde, deltaSigma, lLib] = libraryGeneration(EEm,DEm,ETm,DTm,tm,lf,nLib,Giicm_GG,Sm_GG,gm_GG,Giicm_GC,Sm_GC,gm_GC,Giicm_CC,Sm_CC,gm_CC,lf_SPSL,lf_hybrid)

%% Interpolation parameters
deltaSigma=3; % (MPa)
Npercents=1.5;  % (%)
lflog=logspace(log10(lf/nLib/4),log10(lf),nLib);
lLib=[0 lflog];

%% AA
% FJ - glass-glass overlap
AA{1}=[0 0;100 0];
% FJ - why is this here? for debugging purposes?
i=1;
% FJ - loop through overlap lengths from smallest to largest
% for l=lf/nLib:lf/nLib:lf
for i=1:nLib
    l=lflog(i); %overlap length
    i=i+1;
    
    % Collection of SS curve
    %[SS,S,e,SS_steps]=general_loop_2Stiffness_nonIterative(EEm(1,1),EEm(1,1)/10^6,ETm(1,1),ETm(1,1)/10^6,tm,l);
    % FJ length of process zone may be incorrect (l sent to model) l/4
    % originally
%     Giicm_GG=2;

% [SS_steps,prop,Giicm]=SPSL_VariableGiicSm(EEm(1,1)/2,ETm(1,1)/2,tm,l/2,Giicm_GG,Sm_GG,gm_GG); % Collection of SS curves
[SS_steps,prop,Giicm]=SPSL_VariableGiicSm(EEm(1,1)/2,ETm(1,1)/2,tm,l/2*lf_SPSL,Giicm_GG,Sm_GG,gm_GG); % Collection of SS curves
    AAE(i)=prop(2); % Collection of stiffness
    Epsilon=SS_steps(:,1);
    Sigma=SS_steps(:,2);
    
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
    
    % Interpolation on sigma
    NewSigmaGrid=[0:deltaSigma:smax]';
    
    [NewEpsilon,indUnique]=unique(floor(10^5*NewEpsilon)/10^5);
    NewSigma=NewSigma(indUnique);
    
    NewEpsilonGrid = interp1(NewSigma,NewEpsilon,NewSigmaGrid);
    NewEpsilonGrid(end)=emax;
    
    % Add unloading depending on toughnessclo
    GicTemp=trapz(NewEpsilonGrid,NewSigmaGrid);
    
    % FJ CHECK THIS! Should it be double stiffness?
    dG=Giicm/(ETm(1,1))*100-GicTemp;
    depsilon=2*dG/NewSigmaGrid(end);
    
    % Storage in library
    AA{i}=[NewEpsilonGrid NewSigmaGrid];
    AAde(i)=depsilon;
    
    % Plot - test
    %     figure(19)
    %     hold on
    %     plot(Epsilon,Sigma)
    %     plot(emin,smin,'o')
    %     plot(emax,smax,'o')
    %     plot(NewEpsilon,NewSigma,'--')
    %     plot(NewEpsilonGrid,NewSigmaGrid,'*')
end

% AA{1}=[0 0;100 0];
% i=1;
% % for l=lf/nLib:lf/nLib:lf
% for i=1:nLib
%     l=lflog(i);
%     i=i+1;
%     
%     % Collection of SS curve
%     %%%%%% ?????? [SS,EEE,S,e,SS_steps,Giicm]=general_loop_2Stiffness_nonIterative(EEm(1,2),DEm(1,2),ETm(1,2),DTm(1,2),tm,2*l);
%     % FJ length of process zone may be incorrect (l sent to model) 2*l
%     % originally
% %     Giicm_CC=2;
%     [SS,EEE,S,e,SS_steps,Giicm]=general_loop_2Stiffness_nonIterative_VariableGiicSm(EEm(1,1),DEm(1,1),ETm(1,1),DTm(1,1),tm,2*l,Giicm_CC,Sm_CC,gm_CC);
%     AAE(i)=EEE*1000; % Collection of stiffness
%     Epsilon=SS_steps(:,1);
%     Sigma=SS_steps(:,2);
%     
%     % Correction for plateau
%     % FJ - indBB - is that supposed to be BB or hasn't it changed between
%     % interface types?
%     [M,indBB]=max(Sigma); % Sigma Max
%     indMin=min(find(Sigma>max(Sigma)*(1-Npercents/100)));
%     emin=Epsilon(indMin); % epsilon min in N%
%     smin=max(Sigma)*(1-Npercents/100);
%     smin=Sigma(indMin);
%     indMax=max(find(Sigma>max(Sigma)*(1-Npercents/100)));
%     emax=max(max(Epsilon(indMin:indMax)),Epsilon(indMax)); % epsilon max in N%
%     smax=max(Sigma);
%     
%     % New curve with loading only
%     if smin==smax
%         NewSigma=[Sigma(1:indMin-1) ; smax];
%         NewEpsilon=[Epsilon(1:indMin-1) ; emax];
%     else
%         NewSigma=[Sigma(1:indMin-1) ; smin ; smax];
%         NewEpsilon=[Epsilon(1:indMin-1) ; emin ; emax];
%     end
%     
%     % Interpolation on sigma
%     NewSigmaGrid=[0:deltaSigma:smax]';
%     
%     [NewEpsilon,indUnique]=unique(floor(10^5*NewEpsilon)/10^5);
%     NewSigma=NewSigma(indUnique);
%     
%     if i>123
%         1;
%     end
%     % FJ - this isn't in SPSL code and is giving NaN's. Suspicious else condition is also being satisfied - does this still apply for glass?  - MUST INVESTIGATE!
%     if NewEpsilon<3000;
%         
%         NewEpsilonGrid = interp1(NewSigma,NewEpsilon,NewSigmaGrid);
%         NewEpsilonGrid(end)=emax;
%         
%         % Add unloading depending on toughness
%         GicTemp=trapz(NewEpsilonGrid,NewSigmaGrid);
%         dG=Giicm/(ETm(1,2))*100-GicTemp;
%         % OR dG=Giicm/(ETm(1,2)+tm)*100-GicTemp;
%         depsilon=2*dG/NewSigmaGrid(end);
%         
%         % Storage in library
%         AA{i}=[NewEpsilonGrid NewSigmaGrid];
%         AAde(i)=depsilon;
%         
%     else
%         AA{i}=AA{i-1};
%         AAde(i)=AAde(i-1);
%     end
%     
%     % Plot - test
%     %     figure
%     %     hold on
%     %     plot(Epsilon,Sigma)
%     %     plot(emin,smin,'o')
%     %     plot(emax,smax,'o')
%     %     plot(NewEpsilon,NewSigma,'--')
%     %     plot(NewEpsilonGrid,NewSigmaGrid,'*')
% end
% 
% AA;

%% AB
AB{1}=[0 0;100 0];
i=1;
% for l=lf/nLib:lf/nLib:lf
for i=1:nLib
    l=lflog(i);
    i=i+1;
    
    % Collection of SS curve
    %%%%%% ?????? [SS,EEE,S,e,SS_steps,Giicm]=general_loop_2Stiffness_nonIterative(EEm(1,2),DEm(1,2),ETm(1,2),DTm(1,2),tm,2*l);
    % FJ length of process zone may be incorrect (l sent to model) 2*l
    % originally
%     Giicm_GC=2;
%     [SS,EEE,S,e,SS_steps,Giicm]=general_loop_2Stiffness_nonIterative_VariableGiicSm(EEm(1,2),DEm(1,2),ETm(1,2),DTm(1,2),tm,2*l,Giicm_GC,Sm_GC,gm_GC);
    [SS,EEE,S,e,SS_steps,Giicm]=general_loop_2Stiffness_nonIterative_VariableGiicSm(EEm(1,2),DEm(1,2),ETm(1,2),DTm(1,2),tm,2*l*lf_hybrid,Giicm_GC,Sm_GC,gm_GC);
    ABE(i)=EEE*1000; % Collection of stiffness
    Epsilon=SS_steps(:,1);
    Sigma=SS_steps(:,2);
    
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
    
    % Interpolation on sigma
    NewSigmaGrid=[0:deltaSigma:smax]';
    
    [NewEpsilon,indUnique]=unique(floor(10^5*NewEpsilon)/10^5);
    NewSigma=NewSigma(indUnique);
    
    if i>123
        1;
    end
    % FJ - this isn't in SPSL code and is giving NaN's. Suspicious else condition is also being satisfied - does this still apply for glass?  - MUST INVESTIGATE!
    if NewEpsilon<3000;
         try
             NewEpsilonGrid = interp1(NewSigma,NewEpsilon,NewSigmaGrid);
         catch
             A=NewEpsilon(~isnan(NewEpsilon(1:end-1)))
             B=NewSigma(~isnan(NewSigma(1:end-1)))
             A1=[0];
             B1=[0];
             for j=1:length(B)
                 if max(B(1:j-1))<B(j)
                     B1=[B1 ; B(j)];
                     A1=[A1 ; A(j)];
                 end
             end
             [A1 B1]
             NewSigma=B1;
             NewEpsilon=A1;
             NewSigmaGrid=[0:3:NewSigma(end)]'
         end

    NewEpsilonGrid = interp1(NewSigma,NewEpsilon,NewSigmaGrid);
    NewEpsilonGrid(end)=emax;
    
    
    
    % Add unloading depending on toughness
    GicTemp=trapz(NewEpsilonGrid,NewSigmaGrid);
    dG=Giicm/(ETm(1,2))*100-GicTemp;
    % OR dG=Giicm/(ETm(1,2)+tm)*100-GicTemp;
    depsilon=2*dG/NewSigmaGrid(end);
    
    % Storage in library
    AB{i}=[NewEpsilonGrid NewSigmaGrid];
    ABde(i)=depsilon;
    
     else
         AB{i}=AB{i-1};
         ABde(i)=ABde(i-1);
     end
        % Old libgen
%         NewEpsilonGrid = interp1(NewSigma,NewEpsilon,NewSigmaGrid);
%         NewEpsilonGrid(end)=emax;
%         
%         % Add unloading depending on toughness
%         GicTemp=trapz(NewEpsilonGrid,NewSigmaGrid);
%         dG=Giicm/(ETm(1,2))*100-GicTemp;
%         % OR dG=Giicm/(ETm(1,2)+tm)*100-GicTemp;
%         depsilon=2*dG/NewSigmaGrid(end);
%         
%         % Storage in library
%         AB{i}=[NewEpsilonGrid NewSigmaGrid];
%         ABde(i)=depsilon;
%         
%     else
%         AB{i}=AB{i-1};
%         ABde(i)=ABde(i-1);
%     end
    
    % Plot - test
    %     figure
    %     hold on
    %     plot(Epsilon,Sigma)
    %     plot(emin,smin,'o')
    %     plot(emax,smax,'o')
    %     plot(NewEpsilon,NewSigma,'--')
    %     plot(NewEpsilonGrid,NewSigmaGrid,'*')
end

AB;

%% BB
% carbon-carbon overlap
BB{1}=[0 0;100 0];
i=1;
% for l=lf/nLib:lf/nLib:lf
for i=1:nLib
    l=lflog(i);
    i=i+1;
    
    % Collection of SS curve
    %[SS,S,e,SS_steps]=general_loop_2Stiffness_nonIterative(EEm(2,2),EEm(2,2)/10^6,ETm(2,2),ETm(2,2)/10^6,tm,l);
    % FJ length of process zone may be incorrect (l sent to model)
%     Giicm_CC=2;
%     [SS_steps,prop,Giicm]=SPSL_VariableGiicSm(EEm(2,2)/2,ETm(2,2)/2,tm,l/2,Giicm_CC,Sm_CC,gm_CC);
[SS_steps,prop,Giicm]=SPSL_VariableGiicSm(EEm(2,2)/2,ETm(2,2)/2,tm,l/2*lf_SPSL,Giicm_CC,Sm_CC,gm_CC);
    BBE(i)=prop(2); % Collection of stiffness
    Epsilon=SS_steps(:,1);
    Sigma=SS_steps(:,2);
    
    
    
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
    
    % Interpolation on sigma
    NewSigmaGrid=[0:deltaSigma:smax]';
    
    [NewEpsilon,indUnique]=unique(floor(10^5*NewEpsilon)/10^5);
    NewSigma=NewSigma(indUnique);
    
    NewEpsilonGrid = interp1(NewSigma,NewEpsilon,NewSigmaGrid);
    NewEpsilonGrid(end)=emax;
    
    % Add unloading depending on toughness
    GicTemp=trapz(NewEpsilonGrid,NewSigmaGrid);
    dG=Giicm/(ETm(2,2))*100-GicTemp;
    depsilon=2*dG/NewSigmaGrid(end);
    
    % Storage in library
    BB{i}=[NewEpsilonGrid NewSigmaGrid];
    BBde(i)=depsilon;
    
    % Plot - test
    %     figure
    %     hold on
    %     plot(Epsilon,Sigma)
    %     plot(emin,smin,'o')
    %     plot(emax,smax,'o')
    %     plot(NewEpsilon,NewSigma,'--')
    %     plot(NewEpsilonGrid,NewSigmaGrid,'*')
end

% BB{1}=[0 0;100 0];
% i=1;
% % for l=lf/nLib:lf/nLib:lf
% for i=1:nLib
%     l=lflog(i);
%     i=i+1;
%     
%     % Collection of SS curve
%     %%%%%% ?????? [SS,EEE,S,e,SS_steps,Giicm]=general_loop_2Stiffness_nonIterative(EEm(1,2),DEm(1,2),ETm(1,2),DTm(1,2),tm,2*l);
%     % FJ length of process zone may be incorrect (l sent to model) 2*l
%     % originally
% %     Giicm_CC=2;
%     [SS,EEE,S,e,SS_steps,Giicm]=general_loop_2Stiffness_nonIterative_VariableGiicSm(EEm(2,2),DEm(2,2),ETm(2,2),DTm(2,2),tm,2*l,Giicm_CC,Sm_CC,gm_CC);
%     BBE(i)=EEE*1000; % Collection of stiffness
%     Epsilon=SS_steps(:,1);
%     Sigma=SS_steps(:,2);
%     
%     % Correction for plateau
%     % FJ - indBB - is that supposed to be BB or hasn't it changed between
%     % interface types?
%     [M,indBB]=max(Sigma); % Sigma Max
%     indMin=min(find(Sigma>max(Sigma)*(1-Npercents/100)));
%     emin=Epsilon(indMin); % epsilon min in N%
%     smin=max(Sigma)*(1-Npercents/100);
%     smin=Sigma(indMin);
%     indMax=max(find(Sigma>max(Sigma)*(1-Npercents/100)));
%     emax=max(max(Epsilon(indMin:indMax)),Epsilon(indMax)); % epsilon max in N%
%     smax=max(Sigma);
%     
%     % New curve with loading only
%     if smin==smax
%         NewSigma=[Sigma(1:indMin-1) ; smax];
%         NewEpsilon=[Epsilon(1:indMin-1) ; emax];
%     else
%         NewSigma=[Sigma(1:indMin-1) ; smin ; smax];
%         NewEpsilon=[Epsilon(1:indMin-1) ; emin ; emax];
%     end
%     
%     % Interpolation on sigma
%     NewSigmaGrid=[0:deltaSigma:smax]';
%     
%     [NewEpsilon,indUnique]=unique(floor(10^5*NewEpsilon)/10^5);
%     NewSigma=NewSigma(indUnique);
%     
%     if i>123
%         1;
%     end
%     % FJ - this isn't in SPSL code and is giving NaN's. Suspicious else condition is also being satisfied - does this still apply for glass?  - MUST INVESTIGATE!
%     if NewEpsilon<3000;
%         
%         NewEpsilonGrid = interp1(NewSigma,NewEpsilon,NewSigmaGrid);
%         NewEpsilonGrid(end)=emax;
%         
%         % Add unloading depending on toughness
%         GicTemp=trapz(NewEpsilonGrid,NewSigmaGrid);
%         dG=Giicm/(ETm(1,2))*100-GicTemp;
%         % OR dG=Giicm/(ETm(1,2)+tm)*100-GicTemp;
%         depsilon=2*dG/NewSigmaGrid(end);
%         
%         % Storage in library
%         BB{i}=[NewEpsilonGrid NewSigmaGrid];
%         BBde(i)=depsilon;
%         
%     else
%         BB{i}=BB{i-1};
%         BBde(i)=BBde(i-1);
%     end
%     
%     % Plot - test
%     %     figure
%     %     hold on
%     %     plot(Epsilon,Sigma)
%     %     plot(emin,smin,'o')
%     %     plot(emax,smax,'o')
%     %     plot(NewEpsilon,NewSigma,'--')
%     %     plot(NewEpsilonGrid,NewSigmaGrid,'*')
% end
% 
% BB;

%% Plotting the librabies
% figure
% subplot(1,3,1);
% hold on
% for j=1:length(AA)
%     SS=AA{j};
%     plot(SS(:,1),SS(:,2));
% end
% subplot(1,3,2);
% hold on
% for j=1:length(AB)
%     SS=AB{j};
%     plot(SS(:,1),SS(:,2));
% end
% subplot(1,3,3);
% hold on
% for j=1:length(BB)
%     SS=BB{j};
%     plot(SS(:,1),SS(:,2));
% end

end
