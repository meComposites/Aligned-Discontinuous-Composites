function Pre_defect_stiffness=Find_initial_stiffness(n_rows_outer_layer,n_columns_outer_layer,RVE_outer_layer,E_fibres,MFT1,MFT2,Mt,MSS,GiicmAvg,G,RandomLength,EpsilonGrid,EpsilonMax,Mint,Mfib,Mfail,MNumNeighbour,lf,ds,Vf,tm,inter_fibre_variability)

% Initiation: at the beginning, all interfaces have to be calculated %%%%%%
% FJ - preallocate arrays
ifail=[];
jfail=[];

% FJ - preallocate ifail and jfail vectors with i and j locations of
% failures. First time round you want every fibre to be assumed to be
% failed, so that you evaluate every point. From here on in, only
% recalculate when there is a fibre failure.

% FJ - modified for cracks on outside of RVE
for i=2:n_rows_outer_layer-1
    % for i=1:n_rows
    % ### FJ - na modified to support naxnb rectangle RVE ###
    for j=2:n_columns_outer_layer-1
        %     for j=1:n_columns
        ifail=[ifail, i];
        jfail=[jfail, j];
    end
end



% Interfaces derivation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


iv=[];
jv=[];
for kij=1:length(ifail)
    ifromfail=[2*ifail(kij)-2,2*ifail(kij)-1,2*ifail(kij),2*ifail(kij)-1];
    jfromfail=[jfail(kij),jfail(kij),jfail(kij),jfail(kij)-1];
    for kijff=1:4
        iv=[iv ifromfail(kijff)];
        jv=[jv jfromfail(kijff)];
    end
end

listInt=unique([iv' jv'],'rows');
for iList=1:size(listInt,1)
    
    i=listInt(iList,1);
    j=listInt(iList,2);
    
    
    if1=(i+(mod(i,2)==1))/2;
    jf1=j;
    if2=(i+1+(mod(i,2)==0))/2;
    jf2=(j+(mod(i,2)==1));
    
    
    if if1<=0 || jf1<=0 || if2<=0 || jf2<=0 || mod(if1,1)~=0 || mod(jf1,1)~=0 || mod(if2,1)~=0 || mod(jf2,1)~=0
        1;
    end
    
    type1=RVE_outer_layer(if1,jf1);
    type2=RVE_outer_layer(if2,jf2);
    
    E1=E_fibres(if1,jf1)*1.00001;
    E2=E_fibres(if2,jf2);
    EE_int=E1+E2;
    DE_int=E2-E1;
    
    T1=MFT1(i,j)*1.00001;
    T2=MFT2(i,j);
    ET_int=T1+T2;
    DT_int=T2-T1;
    
    if E2*T2==E1*T1
        signDeltaK=1;
    else
        signDeltaK=abs(E2*T2-E1*T1)/(E2*T2-E1*T1);
    end
    
    switch inter_fibre_variability
        case 1
            tm_int=Mt(i,j);
        otherwise
            tm_int=tm;
    end
    
    Sm_int=MSS(i,j);
    
    Icracks=findIcrack(if1,jf1,if2,jf2,RandomLength,Mfail,lf);
    
    nCrack=size(Icracks,2)-1;
    clearvars SSint
    
    
    
    for k=1:nCrack
        lk=Icracks(1,k+1)-Icracks(1,k);
        
        [xInt,yInt,EInt,deInt]=SSInterface_EffectsOfRandomness(EE_int,DE_int,ET_int,DT_int,signDeltaK,tm_int,Sm_int,GiicmAvg,G,lk,lf);
        
        %                 xInt=xInt/3;
        SSint{k}={xInt,yInt,EInt,deInt,lk};
        
    end
    
    [xSeries,ySeries]=CombineInSeries(SSint,ds,lf,EpsilonMax); % OKKKKKKK!
    
    %         if min(ySeries(2:end-1))==0||max(xSeries>1E+03)
    %             1;
    %         end
    
    [sigmaij] = EpsilonInterpolation(xSeries,ySeries,EpsilonGrid,EpsilonMax);
    Mint(i,j,:)=sigmaij;
    
end


Mint([2 end-1],[2:end-1],:)=0; % removing ints on the edges
Mint([2:end-1],1,:)=0;
Mint([1:2:end-1],end-1,:)=0;

% Fibre derivation for overal response %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Good
Fibre_stress=bsxfun(@rdivide,(Mint(2:2:(end-2),2:(end-1),2)+Mint(3:2:(end-1),2:(end-1),2)+Mint(4:2:end,2:(end-1),2)+Mint(3:2:(end-1),1:(end-2),2)),MNumNeighbour(2:end-1,2:end-1));
Pre_defect_stiffness=Fibre_stress/EpsilonGrid(2)*Vf;
end
%% Epsilon interpolation
    function [SigmaIntOnEpsilon] = EpsilonInterpolation(x,y,EpsilonGrid,EpsilonMax)
        
%         if min(y(2:end-1))==0||max(x>1E+03)
%             1;
%         end
        
        SigmaIntOnEpsilon = interp1([x ; EpsilonMax+10^(-6)],[y ; 0],EpsilonGrid);
        %         SigmaIntOnEpsilon = interp1([x ; x(end)+10^(-6)],[y ; 0],EpsilonGrid);
        
    end

%% Find list/position of cracks in an interface
    function Icracks=findIcrack(if1,jf1,if2,jf2,RandomLength,Mfail,lf)
        
        Icracks1=[RandomLength(if1,jf1) Mfail{if1,jf1}];
        Icracks1(2,1:end)=1;
        Icracks2=[RandomLength(if2,jf2)+lf*(RandomLength(if2,jf2)<RandomLength(if1,jf1))...             % Natural edge of the fibre
            Mfail{if2,jf2}+lf*(Mfail{if2,jf2}<RandomLength(if1,jf1))-lf*(Mfail{if2,jf2}>RandomLength(if1,jf1)+lf)];   % if left/right of fibre 1
        Icracks2(2,1:end)=2;
        
        Icracks=[Icracks1 Icracks2];
        Icracks=sortrows(Icracks',1)';
        
        Icracks=[Icracks [RandomLength(if1,jf1)+lf;1]];
        
    end

%% Collection of interface properties
%     function [xInt,yInt,EInt,deInt]=SSInterface(type1,type2,ndl)
    function [xInt,yInt,EInt,deInt]=SSInterface(type1,type2,length)
        
        if libGen == 2
            %% %%%%% FJ - SORAIA: this section of the fucntion will need to change when random diameters are added
            % %% for when calculating interfaces individually
            if type1+type2==0 % Glass- Glass
                [xInt,yInt,EInt,deInt] = SS_HBaM(EEm(1,1),abs(DEm(1,1)),ETm(1,1),abs(DTm(1,1)),tm,length,SmAvg,GiicmAvg,G);
            elseif type1+type2==1 % Glass - Carbon
                [xInt,yInt,EInt,deInt] = SS_HBaM(EEm(1,2),abs(DEm(1,2)),ETm(1,2),abs(DTm(1,2)),tm,length,SmAvg,GiicmAvg,G);
            else % Carbon - Carbon
                [xInt,yInt,EInt,deInt] = SS_HBaM(EEm(2,2),abs(DEm(2,2)),ETm(2,2),abs(DTm(2,2)),tm,length,SmAvg,GiicmAvg,G);
            end
        else
            ndl=length;
            % for when calling library %%
            if type1+type2==0 % Glass- Glass
                SSI=AA{ndl};
                EInt=AAE(ndl);  % stiffness of "big length"
                deInt=AAde(ndl); % depsilon of "small length"
            elseif type1+type2==1 % Glass - Carbon
                SSI=AB{ndl};
                EInt=ABE(ndl);  % stiffness of "big length"
                deInt=ABde(ndl); % depsilon of "small length"
            else % Carbon - Carbon
                SSI=BB{ndl};
                EInt=BBE(ndl);  % stiffness of "big length"
                deInt=BBde(ndl); % depsilon of "small length"
            end
            
            xInt=SSI(:,1);
            yInt=SSI(:,2);
        end
        

        
    end
%% Collection of interface properties for Effects of Randomness study
    function [xInt,yInt,EInt,deInt]=SSInterface_EffectsOfRandomness(EE_int,DE_int,ET_int,DT_int,signDeltaK,tm_int,Sm_int,GiicmAvg,G,lk,lf)
        
        if lk <= lf/512
            xInt=[0;100];
            yInt=[0;0];
            EInt=0;
            deInt=0;
        elseif tm_int <= 5E-04
            [EInt,deInt,xInt,yInt] = HBaMBiLinSP_robust_FJPostProcess_update_AvgIntervals(EE_int,signDeltaK*DE_int,ET_int,signDeltaK*DT_int,5E-04,lk,Sm_int,GiicmAvg,G,1000);
        else
            [EInt,deInt,xInt,yInt] = HBaMBiLinSP_robust_FJPostProcess_update_AvgIntervals(EE_int,signDeltaK*DE_int,ET_int,signDeltaK*DT_int,tm_int,lk,Sm_int,GiicmAvg,G,1000);
%         xInt=xInt*100;
%         [xInt,yInt,EInt,deInt] = SS_HBaM(EE_int,signDeltaK*DE_int,ET_int,signDeltaK*DT_int,tm_int,lk,Sm_int,GiicmAvg,G);
        end

    end
%% Combination in series
    function [xSeries,ySeries]=CombineInSeries(SSn,ds,lf,EpsilonMax) %{xInt,yInt,EInt,deInt,lk} in each
        
        n=length(SSn);
        
        Smaxn=[];
        
        for ki=1:n
            Smaxn=[Smaxn max(SSn{ki}{2})];
        end
        % FJ - finds minimum strength of two overlaps
        [Smax,iMin]=min(Smaxn);
        % changed with SP updates to effects of randomness
        %         nSmax=Smax/ds;
        % number of increments equals minimum strength / ds
        nSmax=floor(Smax/ds);
        
        % if less than ds, set to zero
        if nSmax==0
            xSeries=[0 EpsilonMax]';
            ySeries=[0 0]'; % big assumption here ...
        else
            
            % Loading
            
            u=zeros(nSmax+1,1); % Total displacement
            for kl=1:n
                % displacement equals previous displacements + strain times
                % length of overlap
                u=u+SSn{kl}{1}(1:nSmax+1)/100*SSn{kl}{5}; % in mm
            end
            % find epsilon from displacement of total interaction
            Epsilon=u/lf*100; % in (%)
            
            % Unloading
            du=SSn{iMin}{4}/100*SSn{iMin}{5}; % in (mm)
            for ku=[1:iMin-1 iMin+1:n]
                du=du-Smax/(SSn{ku}{3})*SSn{ku}{5}; % in (mm) (and E in MPa)
            end
            dEpsilon=max(10^(-6),du/lf*100); % in (%)
            
            Epsilon(end+1,1)=Epsilon(end,1)+dEpsilon; % in (%)
            
            xSeries=Epsilon;
            ySeries=[0:ds:ds*nSmax 0]';
            
            if min(ySeries(2:end-1)==0)
                1;
            end
            
        end
        
    end

