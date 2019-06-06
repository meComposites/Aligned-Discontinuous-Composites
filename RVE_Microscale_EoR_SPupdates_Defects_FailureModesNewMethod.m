
% function [Ecomposite,MSScurves,alphag,alphac,MminOverlap,MminOverlapType] = HiPerDiF_ModelFullSSStrainStep_RectangleXSection(Vf,lf,pc,pg,HalfLengthC,HalfLengthG,dg,Eg,Sg,dc,Ec,Sc,Vc,G,RVEg,first_fibre_overlap,n_rows,n_columns,HiPerDiF_fuzzyness,Fuzzy_Factor,x_section_arrangement,max_strain,GiicmAvg,SmAvg,gmAvg,GiicmAvg,SmAvg,gmAvg,GiicmAvg,SmAvg,gmAvg,c_width,n_columns_c,mcWei,mgWei,lcWei,lgWei,Sc_frag_const,lf_SPSL,lf_hybrid,AA,AAE,AAde, AB,ABE,ABde, BB,BBE,BBde, ds, lLib)
function [Ecomposite,MSScurves,alphag,alphac,MminOverlap,MminOverlapType,crack_count,RVE_data,tm] = RVE_Microscale_EoR_SPupdates_Defects_FailureModesNewMethod(RVE_loop_counter,deltaEpsilon,Vf,lf,pc,pInterface,pcFailed,pg,pgPartial,pgFailed,HalfLengthC,HalfLengthG,dg,EgAvg,EgCoV,Sg,dc,EcAvg,EcCoV,Sc,Vc,G,RVEg,first_fibre_overlap,n_rows,n_columns,HiPerDiF_fuzzyness,Fuzzy_Factor,x_section_arrangement,max_strain,GiicmAvg,SmAvg,SmCoV,gmAvg,c_width,n_columns_c,mcWei,mgWei,lcWei,lgWei,Sc_frag_const,lf_SPSL,lf_hybrid,libGen,AA,AAE,AAde, AB,ABE,ABde, BB,BBE,BBde, ds, lLib,rng_switch,failure_modes_switch)
% FJ - This function calculates the stress-strain response of a
% representative volume element (RVE) given the following inputs:
% Vf: volume fraction of the composite (all fibre types)
% lf: length of fibres
% dg: diameter of fibre 1 (originally glass fibre)
% Eg: stiffness of fibre 1
% Sg: strength of fibre 1
% dc: diameter of fibre 2 (originally carbon fibre)
% Ec: stiffness of fibre 2
% Sc: strength of fibre 2
% Vc: volume fraction of carbon fibres
% G: shear stiffness of matrix
% RVEg: unknown parameter???
% RandomLengthg: matrix of lengths of fibres at intersection with RVE bundary
inter_fibre_variability = 1;
fibre_modulus_variability = 1;
matrix_strength_variability = 1;
% references:
% [1] - "Prediction of Mechanical Properties of Hybrid Discontinuous
% Composites", J Henry, ECCM17, 2016

% set rng start point (needs to change consistently for every RVE during parallelisation)
switch rng_switch
    case 0
        rng(RVE_loop_counter);
    case 1
        rng('shuffle','simdTwister')       
    otherwise
        rng(rng_switch);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RVE_data=RVE_dataset;
RVE_data.Fibre_stagger_distance=first_fibre_overlap(2:end-1,2:end-1);
% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FJ - calculate equivalent platelet thickness on one side of fibre
% boundary
Tc=dc/4; % equivalent thickness of carbon
Tg=dg/4; % equivalent thickness of glass


t1=[0 0;1 1];
t2=[0 1;0 1];
% FJ - calculate sigma(Ea, Eb) for different  fibre types (c-c,c-g;c-g,g-g)
EEm=EcAvg*(t1+t2)+EgAvg*(2-t1-t2);
% FJ - calculate delta(Ea, Eb) for different  fibre types (c-c,c-g;c-g,g-g)
DEm=(EcAvg-EgAvg)*(t1-t2);
% FJ - calculate sigma(Ta, Tb) for different  fibre types (c-c,c-g;c-g,g-g)
%% %%%%%% FJ - SORAIA: this may change when matrix of diameters is added %%%
ETm=Tc*(t1+t2)+Tg*(2-t1-t2);
% FJ - calculate delta(Ta, Tb) for different  fibre types (c-c,c-g;c-g,g-g)
DTm=(Tc-Tg)*(t1-t2);


%% RVE generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FJ - calculate glass fibre volume fraction (relative to # carbon fibres)
Vg=1-Vc;
% ### FJ - na modified to support naxnb rectangle RVE ###
nf=n_rows.*n_columns;        % Total number of fibres
% FJ - calculate number of carbon fibres, based on area fraction
nc=floor(nf*(Vc/dc^2)/(Vc/dc^2+Vg/dg^2));

% FJ - call RVE_config to generate RVE x-section fibre layout
% ### FJ - na modified to support naxnb rectangle RVE ###
[RVE,nc]=RVE_config_RectangleXSection_HiPerDiF_FailureModes(RVE_loop_counter,n_rows,n_columns,nc,Vc,x_section_arrangement,c_width,n_columns_c,HiPerDiF_fuzzyness,Fuzzy_Factor,rng_switch);
RVE_data.RVE_fibre_type=RVE;
% RVE=[0 0 1 0 0;0 0 1 0 0;0 1 1 1 0;0 0 1 0 0;0 0 1 0 0];
% nc=5;

if abs((nc-(floor(nf*(Vc/dc^2)/(Vc/dc^2+Vg/dg^2))))/(floor(nf*(Vc/dc^2)/(Vc/dc^2+Vg/dg^2))))>0.05
    disp('carbon volume fraction does not match with HiPerDiF arrangement');
end

% FJ - calculate number of glass fibres
ng=nf-nc;

% Generation of the matrix thickness (used for library generation)
%% %%%%%% FJ - SORAIA: this may change when matrix of diameters is added %%%
tm=(sqrt(-Vf*(-((nc/nf)+(ng/nf))*((nc/nf)*dc^2+(ng/nf)*dg^2)*pi+(nc/nf)*Vf*(ng/nf)*(dc-dg)^2*pi)*pi)+(-(nc/nf)*dc-(ng/nf)*dg)*Vf*pi)/(Vf*pi*((nc/nf)+(ng/nf)));

% FJ - add outer layer to implement edge BCs
RVE_outer_layer=[RVE(1,1) RVE(1,:) RVE(1,end);RVE(:,1) RVE RVE(:,end);RVE(end,1) RVE(end,:) RVE(end,end)];
% RVE_outer_layer=ones(18,281);
n_rows_outer_layer=n_rows+2;
n_columns_outer_layer=n_columns+2;

%% Caluculation matrices derivations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FJ - uses the fact that RVE is made from 1s and zeros to generate
% matrices for the hybrid RVE
ME=RVE_outer_layer.*EcAvg+(1-RVE_outer_layer).*EgAvg;         % matrice of average moduli

%% %%%%%%%%%%%%% FJ - SORAIA: this (and the following lines) may need to change to implement different diameters (and spacings)
MD=RVE_outer_layer.*dc+(1-RVE_outer_layer).*dg;         % matrice of diameters

MS=RVE_outer_layer.*Sc+(1-RVE_outer_layer).*Sg;         % matrice of strengths

%Matrix of fibre-moduli, used for each interaction
switch fibre_modulus_variability
    case 1
        E_fibres=defineEf(RVE_outer_layer,EgAvg,EgCoV,EcAvg,EcCoV);
    otherwise
        E_fibres=ME;
end
RVE_data.Fibre_stiffness=E_fibres(2:end-1,2:end-1);

nfc=sum(sum(RVE_outer_layer))/n_rows_outer_layer*n_columns_outer_layer;

%Matrices of interaction thickness
[tv,th]=definetI(n_rows_outer_layer,n_columns_outer_layer,Vf,nfc,dg,dc);
Mt=zeros(2*n_rows_outer_layer-1,n_columns_outer_layer);
Mt(1:2:end,1:end-1)=th;
Mt(2:2:end,1:end)=tv;

RVE_data.Inter_fibre_distance_record=Record_interaction_properties(tv,th,'record');

%Matrices of interaction shear strength


switch matrix_strength_variability
    case 1
        MSS=zeros(2*n_rows_outer_layer-1,n_columns_outer_layer);
        [Sv,Sh]=defineSI(n_rows_outer_layer,n_columns_outer_layer,SmAvg,SmCoV);
        MSS(1:2:end,1:end-1)=Sh;
        MSS(2:2:end,1:end)=Sv;
        RVE_data.Matrix_strength_record=Record_interaction_properties(Sv,Sh,'record');
        
    otherwise
        MSS=ones(2*n_rows_outer_layer-1,n_columns_outer_layer)*SmAvg;
end

switch inter_fibre_variability
    case 1
        %Matrix of fibre thicknesses for interactions with k=1:4 neighbours
        Tf = aux_defineTfI(n_rows_outer_layer,n_columns_outer_layer,MD,tv,th);
    otherwise
        Tf = (RVE_outer_layer==1)*dc/4+(RVE_outer_layer==0)*dg/4;
        Tf = repmat(Tf,1,1,4);
end

% Fibre_diameter=2.*sqrt(sum(Tf.^2,3));
% MD=Fibre_diameter;
% RVE_data.Fibre_diameter = Fibre_diameter(2:end-1,2:end-1);

%Mapping fibre matrices to Interactions matrices
[T1v,T2v,T1h,T2h] = aux_mapping_Tf_to_TI(Tf);

% create matrix of fibre 1 and fibre 2 interactions
MFT1=zeros(2*n_rows_outer_layer-1,n_columns_outer_layer);
MFT1(1:2:end,1:end-1)=T1h;
MFT1(2:2:end,1:end)=T1v;

MFT2=zeros(2*n_rows_outer_layer-1,n_columns_outer_layer);
MFT2(1:2:end,1:end-1)=T2h;
MFT2(2:2:end,1:end)=T2v;


%% % ### FJ - na modified to support naxnb rectangle RVE ###
% FJ - calculate M for the four different overlaps
% % ### FJ - na modified to support naxnb rectangle RVE ###
Mcalc11(2*n_rows_outer_layer-1,n_columns_outer_layer)=0; % Top or Left fibre
Mcalc11(1:2:2*n_rows_outer_layer-1,1:n_columns_outer_layer-1)=(MD(:,1:n_columns_outer_layer-1).^2+MD(:,2:n_columns_outer_layer).^2)./MD(:,1:n_columns_outer_layer-1).^2; % Odd half (horizontal interfaces)
Mcalc11(2:2:2*n_rows_outer_layer-2,1:n_columns_outer_layer)=(MD(1:n_rows_outer_layer-1,:).^2+MD(2:n_rows_outer_layer,:).^2)./MD(1:n_rows_outer_layer-1,:).^2; % Even half (vertical interfaces)
%
% % ### FJ - na modified to support naxnb rectangle RVE ###
Mcalc12(2*n_rows_outer_layer-1,n_columns_outer_layer)=0; % Bottom or Right fibre:
Mcalc12(1:2:2*n_rows_outer_layer-1,1:n_columns_outer_layer-1)=(MD(:,1:n_columns_outer_layer-1).^2+MD(:,2:n_columns_outer_layer).^2)./MD(:,2:n_columns_outer_layer).^2; % Odd half (horizontal interfaces)
Mcalc12(2:2:2*n_rows_outer_layer-2,1:n_columns_outer_layer)=(MD(1:n_rows_outer_layer-1,:).^2+MD(2:n_rows_outer_layer,:).^2)./MD(2:n_rows_outer_layer,:).^2; % Even half (vertical interfaces)


% ### FJ - na modified to support naxnb rectangle RVE ###
Mcalc21(2*n_rows_outer_layer-1,n_columns_outer_layer)=0; % Bottom or Right fibre
Mcalc21(1:2:2*n_rows_outer_layer-1,1:n_columns_outer_layer-1)=1+(MD(:,2:n_columns_outer_layer).^2.*ME(:,2:n_columns_outer_layer))./(MD(:,1:n_columns_outer_layer-1).^2.*ME(:,1:n_columns_outer_layer-1)); % Odd half (horizontal interfaces)
Mcalc21(2:2:2*n_rows_outer_layer-2,1:n_columns_outer_layer)=1+(MD(2:n_rows_outer_layer,:).^2.*ME(2:n_rows_outer_layer,:))./(MD(1:n_rows_outer_layer-1,:).^2.*ME(1:n_rows_outer_layer-1,:)); % Even half (vertical interfaces)

% ### FJ - na modified to support naxnb rectangle RVE ###
Mcalc22(2*n_rows_outer_layer-1,n_columns_outer_layer)=0; % Top or Left fibre
Mcalc22(1:2:2*n_rows_outer_layer-1,1:n_columns_outer_layer-1)=1+(MD(:,1:n_columns_outer_layer-1).^2.*ME(:,1:n_columns_outer_layer-1))./(MD(:,2:n_columns_outer_layer).^2.*ME(:,2:n_columns_outer_layer)); % Odd half (horizontal interfaces)
Mcalc22(2:2:2*n_rows_outer_layer-2,1:n_columns_outer_layer)=1+(MD(1:n_rows_outer_layer-1,:).^2.*ME(1:n_rows_outer_layer-1,:))./(MD(2:n_rows_outer_layer,:).^2.*ME(2:n_rows_outer_layer,:)); % Even half (vertical interfaces)
% ### FJ - na modified to support naxnb rectangle RVE ###
% Mcalc21(2*n_rows_outer_layer-1,n_columns_outer_layer)=0; % Bottom or Right fibre
% Mcalc21(1:2:2*n_rows_outer_layer-1,1:n_columns_outer_layer-1)=1+(MD(:,2:n_columns_outer_layer).^2.*E_fibres(:,2:n_columns_outer_layer))./(MD(:,1:n_columns_outer_layer-1).^2.*E_fibres(:,1:n_columns_outer_layer-1)); % Odd half (horizontal interfaces)
% Mcalc21(2:2:2*n_rows_outer_layer-2,1:n_columns_outer_layer)=1+(MD(2:n_rows_outer_layer,:).^2.*E_fibres(2:n_rows_outer_layer,:))./(MD(1:n_rows_outer_layer-1,:).^2.*E_fibres(1:n_rows_outer_layer-1,:)); % Even half (vertical interfaces)
% 
% % ### FJ - na modified to support naxnb rectangle RVE ###
% Mcalc22(2*n_rows_outer_layer-1,n_columns_outer_layer)=0; % Top or Left fibre
% Mcalc22(1:2:2*n_rows_outer_layer-1,1:n_columns_outer_layer-1)=1+(MD(:,1:n_columns_outer_layer-1).^2.*E_fibres(:,1:n_columns_outer_layer-1))./(MD(:,2:n_columns_outer_layer).^2.*E_fibres(:,2:n_columns_outer_layer)); % Odd half (horizontal interfaces)
% Mcalc22(2:2:2*n_rows_outer_layer-2,1:n_columns_outer_layer)=1+(MD(1:n_rows_outer_layer-1,:).^2.*E_fibres(1:n_rows_outer_layer-1,:))./(MD(2:n_rows_outer_layer,:).^2.*E_fibres(2:n_rows_outer_layer,:)); % Even half (vertical interfaces)

% Weibull strength distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFIED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FJ - not required for effects of randomness study, but left to make
% things run (could remove for optimisation)
if HalfLengthC==1
    lfWeiC=lf/2;
else
    lfWeiC=lf;
end

% CARBON %%%%%%%%%%%%%%%%%%%%%
% FJ - Weibull parameters no longer recorded here, but set in Experiment
% wrapper function.
XcWei=Sc;%3430;

XcWei=Sc*(lcWei/lfWeiC)^(1/mcWei);%4344;


lcWei=lfWeiC;

% FJ - sort stochastic strength smallest to largest
FcWei=sort(rand(n_rows_outer_layer,n_columns_outer_layer,4),3);
MXcWei=XcWei/gamma(1+1/mcWei)*(-(mcWei+1)*lcWei/(lfWeiC/4)*log(1-FcWei)*2).^(1/mcWei);

% Run the model for different overlap length (Carbon)- FJ - now includes
% variable matrix properties between carbon and glass
%% %%%%%%%%% FJ - SORAIA: a change might be required here for the diameters - perhaps just work with average values?
[~,~,~,~,~,fieldx,~,~,fieldXB,~,~]=SPSL_VariableGiicSm(EEm(2,2)/2,ETm(2,2)/2,tm,lfWeiC/4,GiicmAvg,SmAvg,gmAvg);
[~,~,~,~,~,fieldx1,~,~,fieldXB1,~,~]=SPSL_VariableGiicSm(EEm(2,2)/2,ETm(2,2)/2,tm,1*lfWeiC/10,GiicmAvg,SmAvg,gmAvg);
[~,~,~,~,~,fieldx2,~,~,fieldXB2,~,~]=SPSL_VariableGiicSm(EEm(2,2)/2,ETm(2,2)/2,tm,2*lfWeiC/10,GiicmAvg,SmAvg,gmAvg);
[~,~,~,~,~,fieldx3,~,~,fieldXB3,~,~]=SPSL_VariableGiicSm(EEm(2,2)/2,ETm(2,2)/2,tm,3*lfWeiC/10,GiicmAvg,SmAvg,gmAvg);
[~,~,~,~,~,fieldx4,~,~,fieldXB4,~,~]=SPSL_VariableGiicSm(EEm(2,2)/2,ETm(2,2)/2,tm,4*lfWeiC/10,GiicmAvg,SmAvg,gmAvg);

% % Run the model for different overlap length (Carbon)
% [~,~,~,~,~,fieldx,~,~,fieldXB,~,~]=SPSL(EEm(2,2)/2,ETm(2,2)/2,tm,lfWeiC/4);
% [~,~,~,~,~,fieldx1,~,~,fieldXB1,~,~]=SPSL(EEm(2,2)/2,ETm(2,2)/2,tm,1*lfWeiC/10);
% [~,~,~,~,~,fieldx2,~,~,fieldXB2,~,~]=SPSL(EEm(2,2)/2,ETm(2,2)/2,tm,2*lfWeiC/10);
% [~,~,~,~,~,fieldx3,~,~,fieldXB3,~,~]=SPSL(EEm(2,2)/2,ETm(2,2)/2,tm,3*lfWeiC/10);
% [~,~,~,~,~,fieldx4,~,~,fieldXB4,~,~]=SPSL(EEm(2,2)/2,ETm(2,2)/2,tm,4*lfWeiC/10);

n=7;
% real complex shape
[xsp, ia, ~]=unique(round(100*[lfWeiC/10+fieldx(n,:) 3*lfWeiC/10+fieldx(n,:)])/100);
XBsp=[fieldXB(n,:) fieldXB(n,end:-1:1)];
XBsp=XBsp(ia);
xsp(isnan(xsp)) = [];
XBsp(isnan(XBsp)) = [];

[xsp1, ia1, ~]=unique(round(100*[lfWeiC/10+fieldx1(n,:)])/100);
XBsp1=[fieldXB1(n,:)];
XBsp1=XBsp1(ia1);
xsp1(isnan(xsp1)) = [];
XBsp1(isnan(XBsp1)) = [];

[xsp2, ia2, ~]=unique(round(100*[2*lfWeiC/10+fieldx2(n,:)])/100);
XBsp2=[fieldXB2(n,:)];
XBsp2=XBsp2(ia2);
xsp2(isnan(xsp2)) = [];
XBsp2(isnan(XBsp2)) = [];

[xsp3, ia3, ~]=unique(round(100*[3*lfWeiC/10+fieldx3(n,:)])/100);
XBsp3=[fieldXB3(n,:)];
XBsp3=XBsp3(ia3);
xsp3(isnan(xsp3)) = [];
XBsp3(isnan(XBsp3)) = [];

[xsp4, ia4, ~]=unique(round(100*[4*lfWeiC/10+fieldx4(n,:)])/100);
XBsp4=[fieldXB4(n,:)];
XBsp4=XBsp4(ia4);
xsp4(isnan(xsp4)) = [];
XBsp4(isnan(XBsp4)) = [];

% interpolation
x=0:lfWeiC/100:lfWeiC;
XB1=interp1(xsp1,XBsp1,x);
XB1(isnan(XB1)) = [];
XB2=interp1(xsp2,XBsp2,x);
XB2(isnan(XB2)) = [];
XB3=interp1(xsp3,XBsp3,x);
XB3(isnan(XB3)) = [];
XB4=interp1(xsp4,XBsp4,x);
XB4(isnan(XB4)) = [];
hold on
% plot(x,XB1);
% plot(x,XB2);
% plot(x,XB3);
% plot(x,XB4);
XBav=([XB1 fliplr(XB4(1:end-1))] + [XB2 fliplr(XB3(1:end-1))]);
XBav=(XBav+fliplr(XBav))/4;
% figure(25)
% plot(x,XBav);

% bilinear
XBlin=0:max(XBav)/50:max(XBav);
XBlin=[XBlin XBlin(end-1:-1:1)];
% plot(x,XBlin);

% FJ - linear stress field
% XBav=ones(1,length(XBav)).*max(XBav);
% lf=25;
% calculate weibull distributions
% FJ - modified to include 1/lf term
Wf=-(lfWeiC/4)/lcWei*1/lfWeiC*trapz(x,(XBav/(XcWei/gamma(1+1/mcWei))).^mcWei);           % Found numerically (real)
Wflin=-(lfWeiC/4)/lcWei*1/(mcWei+1)*(max(XBlin)/(XcWei/gamma(1+1/mcWei))).^mcWei;   % Found analytically (linear)
alphac=Wflin/Wf; % ratio of survival prob (log formula)
% alphac=1;

% Correction on C strength distribution
MXcWei=XcWei/gamma(1+1/mcWei)*(-(mcWei+1)*lcWei/(lfWeiC/4)*log(1-FcWei)*(alphac)).^(1/mcWei);
% lf=4;

% GLASS %%%%%%%%%%%%%%%%%
% lgWei=25;
% m parametric - glass weibull modulus influences strain to failure
% mgWei=15;
% mgWei=5.7119;
% mgWei=3;

% FJ - not necessary for effects of randomness study
if HalfLengthG==1
    lfWeiG=lf/2;
else
    lfWeiG=lf;
end

XgWei=Sg;%4344;

XgWei=Sg*(lgWei/lfWeiG)^(1/mgWei);%4344;
lgWei=lfWeiG;

% FJ - don't sort smallet to largest - keep strength constant
% FgWei=sort(rand(n_rows_outer_layer,n_columns_outer_layer,4),3);
% FgWei(:,:,2:4)=0.5*repmat(FgWei(:,:,1),1,1,3);
% FJ - sort strengths smallest to largest
FgWei=sort(rand(n_rows_outer_layer,n_columns_outer_layer,4),3);

% FgWei=sort(rand(n_rows_outer_layer,n_columns_outer_layer,4),3);
MXgWei=XgWei/gamma(1+1/mgWei)*(-(mgWei+1)*lgWei/(lfWeiG/4)*log(1-FgWei)*2).^(1/mgWei);

% FJ - altered to enable different matrix proiperties between different
% fibre types
%% %%%%%%%%% FJ - SORAIA: this may need changing, perhaps work with average values for calculating stress field? Or not?
[~,~,~,~,~,fieldx,~,~,fieldXB,~,~]=SPSL_VariableGiicSm(EEm(1,1)/2,ETm(1,1)/2,tm,lfWeiG/4,GiicmAvg,SmAvg,gmAvg);
[~,~,~,~,~,fieldx1,~,~,fieldXB1,~,~]=SPSL_VariableGiicSm(EEm(1,1)/2,ETm(1,1)/2,tm,1*lfWeiG/10,GiicmAvg,SmAvg,gmAvg);
[~,~,~,~,~,fieldx2,~,~,fieldXB2,~,~]=SPSL_VariableGiicSm(EEm(1,1)/2,ETm(1,1)/2,tm,2*lfWeiG/10,GiicmAvg,SmAvg,gmAvg);
[~,~,~,~,~,fieldx3,~,~,fieldXB3,~,~]=SPSL_VariableGiicSm(EEm(1,1)/2,ETm(1,1)/2,tm,3*lfWeiG/10,GiicmAvg,SmAvg,gmAvg);
[~,~,~,~,~,fieldx4,~,~,fieldXB4,~,~]=SPSL_VariableGiicSm(EEm(1,1)/2,ETm(1,1)/2,tm,4*lfWeiG/10,GiicmAvg,SmAvg,gmAvg);

% [~,~,~,~,~,fieldx,~,~,fieldXB,~,~]=SPSL(EEm(1,1)/2,ETm(1,1)/2,tm,lfWeiG/4);
% [~,~,~,~,~,fieldx1,~,~,fieldXB1,~,~]=SPSL(EEm(1,1)/2,ETm(1,1)/2,tm,1*lfWeiG/10);
% [~,~,~,~,~,fieldx2,~,~,fieldXB2,~,~]=SPSL(EEm(1,1)/2,ETm(1,1)/2,tm,2*lfWeiG/10);
% [~,~,~,~,~,fieldx3,~,~,fieldXB3,~,~]=SPSL(EEm(1,1)/2,ETm(1,1)/2,tm,3*lfWeiG/10);
% [~,~,~,~,~,fieldx4,~,~,fieldXB4,~,~]=SPSL(EEm(1,1)/2,ETm(1,1)/2,tm,4*lfWeiG/10);

n=7;
% real complex shape
[xsp, ia, ~]=unique(round(100*[lfWeiG/10+fieldx(n,:) 3*lfWeiG/10+fieldx(n,:)])/100);
XBsp=[fieldXB(n,:) fieldXB(n,end:-1:1)];
XBsp=XBsp(ia);
xsp(isnan(xsp)) = [];
XBsp(isnan(XBsp)) = [];

[xsp1, ia1, ~]=unique(round(100*[lfWeiG/10+fieldx1(n,:)])/100);
XBsp1=[fieldXB1(n,:)];
XBsp1=XBsp1(ia1);
xsp1(isnan(xsp1)) = [];
XBsp1(isnan(XBsp1)) = [];

[xsp2, ia2, ~]=unique(round(100*[2*lfWeiG/10+fieldx2(n,:)])/100);
XBsp2=[fieldXB2(n,:)];
XBsp2=XBsp2(ia2);
xsp2(isnan(xsp2)) = [];
XBsp2(isnan(XBsp2)) = [];

[xsp3, ia3, ~]=unique(round(100*[3*lfWeiG/10+fieldx3(n,:)])/100);
XBsp3=[fieldXB3(n,:)];
XBsp3=XBsp3(ia3);
xsp3(isnan(xsp3)) = [];
XBsp3(isnan(XBsp3)) = [];

[xsp4, ia4, ~]=unique(round(100*[4*lfWeiG/10+fieldx4(n,:)])/100);
XBsp4=[fieldXB4(n,:)];
XBsp4=XBsp4(ia4);
xsp4(isnan(xsp4)) = [];
XBsp4(isnan(XBsp4)) = [];

% interpolation
x=0:lfWeiG/100:lfWeiG;
XB1=interp1(xsp1,XBsp1,x);
XB1(isnan(XB1)) = [];
XB2=interp1(xsp2,XBsp2,x);
XB2(isnan(XB2)) = [];
XB3=interp1(xsp3,XBsp3,x);
XB3(isnan(XB3)) = [];
XB4=interp1(xsp4,XBsp4,x);
XB4(isnan(XB4)) = [];
hold on
% plot(x,XB1);
% plot(x,XB2);
% plot(x,XB3);
% plot(x,XB4);
XBav=([XB1 fliplr(XB4(1:end-1))] + [XB2 fliplr(XB3(1:end-1))]);
XBav=(XBav+fliplr(XBav))/4;
% figure(25)
% plot(x,XBav);

% bilinear
XBlin=0:max(XBav)/50:max(XBav);
XBlin=[XBlin XBlin(end-1:-1:1)];
% plot(x,XBlin);

%XBav=ones(1,length(XBav)).*max(XBav);
% lfWeiG=25;
% calculate weibull distributions
Wf=1/lfWeiG*trapz(x,(XBav/(XgWei/gamma(1+1/mgWei))).^mgWei);
Wflin=1/(mgWei+1)*(max(XBlin)/(XgWei/gamma(1+1/mgWei))).^mgWei;
alphag=Wflin/Wf;
% alphag=1;

% Correction on G strength distribution
MXgWei=XgWei/gamma(1+1/mgWei)*(-(mgWei+1)*lgWei/(lfWeiG/4)*log(1-FgWei)*(alphag)).^(1/mgWei);
% lf=4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MODIFIED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FJ - Apply the stochastic strength of the correct fibre type to each
% fibre in the RVE. This includes the stengths of the fibres when boken as well
% MstocSWei=bsxfun(@times,MXcWei,RVE_outer_layer) + bsxfun(@times,MXgWei,(1-RVE_outer_layer));
MstocSWei=bsxfun(@times,MXcWei,RVE_outer_layer) + bsxfun(@times,MXgWei,(1-RVE_outer_layer));
% MstocSWei=repmat(MstocSWei,1,1,4);
RVE_data.MstocS_array=MstocSWei(2:end-1,2:end-1,:);
% FJ - Use the stochastic strengths for the current state of fibres only.
% If no fibers are broken, the first 60x60 array of MstocSWei is used. Once
% a fibre is broken, this part of the array cannto be used. Therefore use
% MstocS with a substitutued value from the 2nd array in the third
% dimension of MstocSWei and continue.
MstocS=MstocSWei(:,:,1);


% libraries derivation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FJ - define global variables for stress-strain curve (XX), stiffness
% (XXE) and strain to failure(***why negative values in library provided?***) (XXde) for each type of material interface. ds = stress
% increment and ILib = overlap length.
% global AA AAE AAde AB ABE ABde BB BBE BBde ds lLib

%% FJ - section commented as libGen now carried out in Experiment wrapper function
%%%%%%%%%%%%%%%%%%%%
% Create libraries %
%%%%%%%%%%%%%%%%%%%%

% %FJ - specify the number of different overlap lengths the library will
% %generate stress-strain curves for
% nLib=128; % number of different lengths available
%  % FJ - call libraryGeneration function to generate libraary data
% [AA,AAE,AAde, AB,ABE,ABde, BB,BBE,BBde, ds, lLib]=libraryGeneration(EEm,DEm,ETm,DTm,tm,lf,nLib,GiicmAvg,SmAvg,gmAvg,GiicmAvg,SmAvg,gmAvg,GiicmAvg,SmAvg,gmAvg,lf_SPSL,lf_hybrid);
%  % FJ - save the data as variables in a .mat file
% save lib-HMC-HSC-G0_8S90L3logn128messingaboutwithJoel_FJ.mat AA AAE AAde AB ABE ABde BB BBE BBde ds lLib

% --- OR ---

%%%%%%%%%%%%%%%%%%
% Load libraries %
%%%%%%%%%%%%%%%%%%

% LIN

% load('lib-Jalal-SG-T30.mat')
% load('lib-HMC-HSC-G.5.mat')
% load('lib-HMC-HSC-G2.mat')
% load('lib-HMC-HSC-G.2S60.mat')
% load('lib-HMC-HSC-G.5S30.mat')
% load('lib-HMC-HSC-G2S30.mat')
% load('lib-HMC-HSC-G.5S80.mat')
% load('lib-HMC-HSC-G1S80.mat')
% load('lib-HMC-HSC-G2S80.mat') %%%OK
% load('lib-HMC-HSC-G4S40.mat')
% load('lib-HMC-HSC-G1S100.mat')
% load('lib-HMC-HSC-G2S100.mat') %%% DIFFERENT %%%
% load('lib-HMC-HSC-G2S100L3log.mat')
% load('lib-HMC-HSC-G4S100L3log.mat')
% load('lib-HMC-HSC-G2S100L6log.mat')
% load('lib-HMC-HSC-G2S100L6.mat')
% load('lib-HMC-HSC-G2S100L6n64.mat')
% load('lib-HMC-HSC.mat')
% load('lib-HMC-HSC-0.1.mat')


% LOG
% load('lib-HMC-HSC-G1S50L3logn128.mat')
% load('lib-HMC-HSC-G05S100L3logn128.mat')
% load('lib-HMC-HSC-G2S150L3logn128.mat')
% load('lib-HMC-HSC-G2S100L3logn128.mat')
% load('lib-HMC-HSC-G2S100L3logn128BIS.mat')
% load('lib-HMC-HSC-G2S100L3logn128TER.mat') %%% Joel's BEST %%%
% load('lib-HMC-HSC-G2S100L3logn128TER_FJ.mat')
% load('lib-HMC-HSC-G2S50L3logn128.mat')
% load('lib-HMC-HSC-G4S100L3logn128.mat')
% load('lib-HMC-HSC-G1S85L3logn128TER_FJ.mat') % FJ - model debugging
% load('lib-HMC-HSC-G1S85L3logn128With2lwithFinleyNextToMe.mat') % FJ - debugging
% load('lib-HMC-HSC-G1S85L3logn128TER_FJ_libgencheck.mat') % FJ - debugging

% load('lib-HMC-HSC-G0_8S90L3logn128messingaboutwithJoel_FJ.mat')
% FJ - Random length of overlap at RVE intersection taken as an input from
% wrapper function, because subsequent RVEs need to have same overlaps at
% fibre locations (for same fibre lengths)
% RandomLength=rand(na)*lf;
RandomLength=first_fibre_overlap;

%FJ - plt the carbon / glass fibres and plot the overlap lengths for each
%fibre
% plotRVE()

% Strain step definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FJ - constant strain analysis

di=[-1 0 1 0];
dj=[0 1 0 -1];
% FJ - define strain increment for each step - consider changing?
% deltaEpsilon=0.02;  % (%)
% FJ - define maximum strain [%] (may need to adjust this for different
% materials)
% FJ - MAKRE SURE THIS IS MUCH BIGGER THAN STRAIN TO FAILURe OF THE HIGHEST FAILURE STRAIN MATERIAL
EpsilonMax=max_strain; % (%)
% FJ - number of strain increments = max strain / strain increment
nEpsilon=floor(EpsilonMax/deltaEpsilon)+1;
% FJ - define the strain scale from 0 to max strain, with equal increments
EpsilonGrid=[0:deltaEpsilon:EpsilonMax]';
% FJ - Set failure indicator to 1. Used later to control iterations ***
% return to this ***
indFailure=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFIED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrices initiation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ### FJ - na modified to support naxnb rectangle RVE ###
Mint(2*n_rows_outer_layer-1,n_columns_outer_layer,nEpsilon)=0;
Mint_record=zeros(2*n_rows_outer_layer-1,n_columns_outer_layer,nEpsilon);
% ### FJ - na modified to support naxnb rectangle RVE ###
Mfib(n_rows_outer_layer,n_columns_outer_layer,nEpsilon)=0;
MFailureOverall=Mfib;
% ### FJ - na modified to support naxnb rectangle RVE ###
Mfail{n_rows_outer_layer,n_columns_outer_layer}=[];
MSScurves=Mfib;
% FJ - at the start of the simulation, all fibres are assumed to be
% unbroken
% ### FJ - na modified to support naxnb rectangle RVE ###
MNumNeighbour=4*ones(n_rows_outer_layer,n_columns_outer_layer);

% remove neighbours on edges
MNumNeighbour([2:end-1],[2 end-1])=3;
MNumNeighbour([2 end-1],[2:end-1])=3;
MNumNeighbour([2 end-1],[2 end-1])=2;
MNumNeighbourIni=MNumNeighbour;


% add p (%) of crack in carbon (optional)

crack_count=[0;0;cell(1,1)];
RVE_data.Carbon_crack_count=[0;0;cell(1,1)];
RVE_data.Glass_crack_count=[0;0;cell(1,1)];
RVE_data.Fibre_crack_count=zeros(n_rows+2,n_columns+2);
RVE_data.Fibre_precrack_count=zeros(n_rows+2,n_columns+2);
RVE_data.Prefailed_interface_count=zeros(2*n_rows_outer_layer-1,n_columns_outer_layer);
RVE_data.Fibre_partialfailed_count=zeros(n_rows+2,n_columns+2);
RVE_data.Fibre_prefailed_count=zeros(n_rows+2,n_columns+2);
if failure_modes_switch==1
    RVE_data.Fragmentation_record=cell(n_rows_outer_layer,n_columns_outer_layer);
end
for i=1:n_rows_outer_layer
    for j=1:n_columns_outer_layer
        RVE_data.Fragmentation_record{i,j}=[1,Inf;2,Inf;3,Inf;4,Inf];
    end
end
OverlapData();
% should be able to switch for different libGen
RVE_data.Pre_defect_stiffness=Find_initial_stiffness(n_rows_outer_layer,n_columns_outer_layer,RVE_outer_layer,E_fibres,MFT1,MFT2,Mt,MSS,GiicmAvg,G,RandomLength,EpsilonGrid,EpsilonMax,Mint,Mfib,Mfail,MNumNeighbour,lf,ds,Vf,tm,inter_fibre_variability);
% data storage for debugging of non-robust HBaM calls
RVE_data.Interface_assert_flag_record=zeros(size(RVE_data.Prefailed_interface_count));
RVE_data.Interface_assert_data_record=cell(size(RVE_data.Prefailed_interface_count));
RVE_data.Fibre_assert_flag_record=zeros(n_rows_outer_layer,n_columns_outer_layer);



%% add defects
ResidualCrackC(pc);
ResidualCrackG(pg);

PreFailedInterface(pInterface)

PreFailedFibres(pcFailed,pgFailed,'Full');
%%
% PreFailedFibres(pcPartial,pgPartial,'Partial');

% %% FJ - not used in effects of randomness study
% % If HiPerDiF half length activated for carbon
% if HalfLengthC==1
%     ResidualCrackC(1);
% end
%
% % If HiPerDiF half length activated for glass
% if HalfLengthG==1
%     ResidualCrackG(1);
% end
%
% % If HiPerDiF half length not activated
% if HalfLengthC==0 && HalfLengthG==0
%     ResidualCrackC(pc);
%     ResidualCrackG(pg);
% end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MODIFIED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FJ - keep while loop until failure. Loop will stop if max number of
% strain increments is reached.
while indFailure < nEpsilon
    indFailurePrevious=indFailure;
    
    % Interfaces derivation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % FJ - previously interfaces were derived for every strain increment,
    % now only evaluated when there is a fibre failure? *** check this
    % All interfaces derived: ----------------------
    % for i=1:2*na-1
    %    % ### na should be modified? ###
    %    for j=1:na-(mod(i,2)==1)
    % ----------------------------------------------
    % OR
    % Only modified interfaces derived -------------
    
    %     %FJ - added to evaluate all interfaces each time
    %     for i=2:n_rows_outer_layer-1
    % % for i=1:n_rows
    %     % ### FJ - na modified to support naxnb rectangle RVE ###
    %         for j=2:n_columns_outer_layer-1
    % %     for j=1:n_columns
    %         ifail=[ifail, i];
    %         jfail=[jfail, j];
    %     end
    % end
    
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
    %     % FJ - modified to enable cracks on edges
    %     listInt=zeros(n_rows*n_columns,2);
    %     for ii=1:n_rows
    %         for jj=1:n_columns
    %             listInt(((ii-1)*n_columns+jj),1)=ii;
    %             listInt(((ii-1)*n_columns+jj),2)=jj;
    %         end
    %     end
    
    listInt=unique([iv' jv'],'rows');
    for iList=1:size(listInt,1)
        
        i=listInt(iList,1);
        j=listInt(iList,2);
        % ----------------------------------------------
        
        if1=(i+(mod(i,2)==1))/2;
        jf1=j;
        if2=(i+1+(mod(i,2)==0))/2;
        jf2=(j+(mod(i,2)==1));
        %         if1=i;
        %         if2=i;
        %         jf1=j;
        %         jf2=j;
        
        if if1<=0 || jf1<=0 || if2<=0 || jf2<=0 || mod(if1,1)~=0 || mod(jf1,1)~=0 || mod(if2,1)~=0 || mod(jf2,1)~=0
            1;
        end
        
        if if1==44 && jf1 == 34
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
%         T1=MFT2(i,j);
%         T2=MFT1(i,j);
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
        
        Icracks=findIcrack(if1,jf1,if2,jf2);
        
        nCrack=size(Icracks,2)-1;
        clearvars SSint
        
        
        
        for k=1:nCrack
            if Icracks(2,k)==Icracks(2,k+1)
                if indFailure ==1
                    indFailure;
                end
                lk=Icracks(1,k+1)-Icracks(1,k);
                deInt=0;
                %% %%%%%%% FJ - SORAIA: these bits may need changing
                %                 Lpz=FloatingLPZ_VariableGiicSm(EEm(type1+1,type2+1),DEm(type1+1,type2+1),ETm(type1+1,type2+1),DTm(type1+1,type2+1),tm,lk/2,GiicmAvg,SmAvg,gmAvg);
                
                EInt=((Icracks(2,k)==1)*E2*(T2*4)^2+(Icracks(2,k)==2)*E1*(T1*4)^2)/((T1*4)^2+(T2*4)^2); % in MPa
                %                 % modified stiffness calcs to investigate effect of
                %                 % floating fibres -
                %                 % FJ - fragmented interface stiffness, based on
                %                 % {Alessi2017} - works for initial stiffness! Need to do the whole thing!
                %                 Eh=(Icracks(2,k)==1)*EEm(type2+1,type2+1)/2+(Icracks(2,k)==2)*EEm(type1+1,type1+1)/2;
                %                 Ef=(Icracks(2,k)==1)*EEm(type1+1,type1+1)/2+(Icracks(2,k)==2)*EEm(type2+1,type2+1)/2;
                %                 th=(Icracks(2,k)==1)*ETm(type2+1,type2+1)/2+(Icracks(2,k)==2)*ETm(type1+1,type1+1)/2;
                %                 tf=(Icracks(2,k)==1)*ETm(type1+1,type1+1)/2+(Icracks(2,k)==2)*ETm(type2+1,type2+1)/2;
                %                 F=(G*(lk/2)^2)/(Eh*th);
                %                 K=(Ef*tf)/(Eh*th);
                %                 alpha_frag=(F*(K+1))/(K);
                %                 beta_frag=F/K;
                %                 E_frag=(alpha_frag*sqrt(alpha_frag))/(sqrt(alpha_frag)*beta_frag+(alpha_frag-beta_frag)*tanh(sqrt(alpha_frag)))*(1/1+tm);
                %                 EInt2=E_frag*Eh;
                %                 if E_frag<10 && E_frag>=1
                %                     Eint=E_frag*Eh;
                %                 else
                %                     EInt=((Icracks(2,k)==1)*EEm(type2+1,type2+1)/2*(ETm(type2+1,type2+1)*2)^2+(Icracks(2,k)==2)*EEm(type1+1,type1+1)/2*(ETm(type1+1,type1+1)*2)^2)/((ETm(type1+1,type1+1)*2)^2+(ETm(type2+1,type2+1)*2)^2); % in MPa
                %                 end
                %                 EInt=((Icracks(2,k)==1)*EEm(type2+1,type2+1)/2*(ETm(type2+1,type2+1)*2)^2+(Icracks(2,k)==2)*EEm(type1+1,type1+1)/2*(ETm(type1+1,type1+1)*2)^2)/((ETm(type1+1,type1+1)*2)^2+(ETm(type2+1,type2+1)*2)^2); % in MPa
                
                %                                 EInt=((Icracks(2,k)==1)*EEm(type2+1,type2+1)/2*(ETm(type2+1,type2+1)*2)^2+(Icracks(2,k)==2)*((lk>(2*Lpz))*((((lk-2*Lpz)/lk)*(EEm(type2+1,type2+1)/2*(ETm(type2+1,type2+1)*2)^2))+((2*Lpz/lk)*EEm(type1+1,type1+1)/2*(ETm(type1+1,type1+1)*2)^2)/((ETm(type1+1,type1+1)*2)^2+(ETm(type2+1,type2+1)*2)^2))+(lk<=(2*Lpz))*((EEm(type1+1,type1+1)/2*(ETm(type1+1,type1+1)*2)^2)/((ETm(type1+1,type1+1)*2)^2+(ETm(type2+1,type2+1)*2)^2)))); % in MPa
                SInt=(Icracks(2,k)==1)*MstocS(if2,jf2)+(Icracks(2,k)==2)*MstocS(if1,jf1);
                %                 SInt=(Icracks(2,k)==1)*MstocS(if2,jf2)+(Icracks(2,k)==2)*((lk>(2*Lpz))*((((lk-2*Lpz)/lk)*(MstocS(if1,jf1)))+((2*Lpz/lk)*MstocS(if2,jf2)))+(lk<=(2*Lpz))*MstocS(if2,jf2));
                
                yInt = reshape([0:ds:1*SInt],[],1);
                xInt=yInt./EInt*100;
                
                SSint{k}={xInt,yInt,EInt,deInt,lk};
            else
                lk=Icracks(1,k+1)-Icracks(1,k);
                
                % FJ - this is not needed unless full analytical fracture
                % toughness is calculated - commented out for speed
                %                 if indFailure==1
                %                     Moverlaps(i,j,k)=lk;
                %                 end
                
                %                 nlk=floor(nlk/0.2);
                %                 if nlk==0
                %                     nlk=1;
                %                 elseif nlk>128
                %                     nlk=128;
                %                 end
                
                %                 if type1==type2
                %                     nlk=floor(nlk/0.2);
                %                     if nlk==0
                %                         nlk=1;
                %                     elseif nlk>128
                %                         nlk=128;
                %                     end
                %                 end
                
                %                 lk=lk/1;
                %
                %                 if type1==type2
                %                     lk=lk/1;
                %                 end
                %
                
                if pInterface>0
                    if RVE_data.Prefailed_interface_count(i,j) == 1
                        lk=0;
                    end
                end
                
                if lk<0 || lk>=lf
                    1;
                end
                if libGen == 2
                    [xInt,yInt,EInt,deInt,assert_flag]=SSInterface_EffectsOfRandomness(EE_int,DE_int,ET_int,DT_int,signDeltaK,tm_int,Sm_int,lk);
                    if assert_flag==1
                        RVE_data.Interface_assert_flag_record(i,j)=RVE_data.Interface_assert_flag_record(i,j)+1;
                        RVE_data.Interface_assert_data_record{i,j}=[RVE_data.Interface_assert_data_record{i,j};[EE_int,DE_int,ET_int,DT_int,signDeltaK,tm_int,Sm_int,lk]];
                        RVE_data.Fibre_assert_flag_record(if1,jf1)=RVE_data.Fibre_assert_flag_record(if1,jf1)+1;
                        RVE_data.Fibre_assert_flag_record(if2,jf2)=RVE_data.Fibre_assert_flag_record(if2,jf2)+1;
                    end
                else
                    % FJ - lton no longer needed (direct calculation of interface
                    % values)
                    nlk=lton(lk);
                    [xInt,yInt,EInt,deInt]=SSInterface(type1,type2,nlk);
                end
                
                %                 xInt=xInt/3;
                SSint{k}={xInt,yInt,EInt,deInt,lk};
                
            end
        end
        
        [xSeries,ySeries]=CombineInSeries(SSint); % OKKKKKKK!
        
        %         if min(ySeries(2:end-1))==0||max(xSeries>1E+03)
        %             1;
        %         end
        
        [sigmaij] = EpsilonInterpolation(xSeries,ySeries);
        Mint(i,j,:)=sigmaij;
        
    end
    % if all fibres each time: ------
    % end
    % -----------------------------------
    
    Mint([2 end-1],[2:end-1],:)=0; % removing ints on the edges
    Mint([2:end-1],1,:)=0;
    Mint([1:2:end-1],end-1,:)=0;
    
    % Fibre derivation for overal response %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Good
    Mfib(2:(end-1),2:(end-1),:)=bsxfun(@rdivide,(Mint(2:2:(end-2),2:(end-1),:)+Mint(3:2:(end-1),2:(end-1),:)+Mint(4:2:end,2:(end-1),:)+Mint(3:2:(end-1),1:(end-2),:)),MNumNeighbour(2:end-1,2:end-1));
%     Mfib(2:(end-1),2:(end-1),:)=bsxfun(@rdivide,(Mint(2:2:(end-2),2:(end-1),:).*Tf(2:end-1,2:end-1,1)+Mint(3:2:(end-1),2:(end-1),:).*Tf(2:end-1,2:end-1,2)+Mint(4:2:end,2:(end-1),:).*Tf(2:end-1,2:end-1,3)+Mint(3:2:(end-1),1:(end-2),:).*Tf(2:end-1,2:end-1,4)),((Mint(2:2:(end-2),2:(end-1),:)>0).*Tf(2:end-1,2:end-1,1)+(Mint(3:2:(end-1),2:(end-1),:)>0).*Tf(2:end-1,2:end-1,2)+(Mint(4:2:end,2:(end-1),:)>0).*Tf(2:end-1,2:end-1,3)+(Mint(3:2:(end-1),1:(end-2),:)>0).*Tf(2:end-1,2:end-1,4)));
%     Mfib(:,:,1)=0;
    
    % NaN from divide by zero if fibre has no neighbours. set Mfib to zero
    % for this fibre, as it cannot transfer stress when no neighbours and
    % if it ever did have neighbours, that will be recorded in MSScurves
    [rowsNoNeighbours,columnsNoNeighbours]=ind2sub([n_rows_outer_layer,n_columns_outer_layer],find((MNumNeighbour==0)));
    Mfib(rowsNoNeighbours,columnsNoNeighbours,:)=0;
    % Wrong
    %Mfib(2:(end-1),2:(end-1),:)=(Mint(2:2:(end-2),2:(end-1),:)+Mint(3:2:(end-1),2:(end-1),:)+Mint(4:2:end,2:(end-1),:)+Mint(3:2:(end-1),1:(end-2),:))./4;
    % Dead fibre
    % Mfib(2:(end-1),2:(end-1),:)=bsxfun(@rdivide,(Mint(2:2:(end-2),2:(end-1),:)+Mint(3:2:(end-1),2:(end-1),:)+Mint(4:2:end,2:(end-1),:)+Mint(3:2:(end-1),1:(end-2),:)),MNumNeighbourIni(2:end-1,2:end-1));
    % Mfib(2:(end-1),2:(end-1),:)=bsxfun(@times,Mfib(2:(end-1),2:(end-1),:),MNumNeighbourIni(2:end-1,2:end-1)==MNumNeighbour(2:end-1,2:end-1));
    
    
    %plotstress()
    
    
    % Derivation of stress peaks in fibre (for failure checking) %%%%%%%%%%
    [Mpeakmax,indMpeakmax]=findMpeak();
    
    % Check failure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %MFailure=repmat(MstocS,1,1,size(Mint,3))<Mpeakmax(:,:,:);
    MFailure=bsxfun(@lt,MstocS,Mpeakmax);
    
    MFailure([1 end],:,:)=0;        % No failure on the edges
    MFailure(:,[1 end],:)=0;        % No failure on the edges
    
    nbrFailure(1:size(Mint,3))=(sum(sum(MFailure)));
    if isempty(find(nbrFailure>0,1,'first'))
        indFailure=indFailure+1; % ????????????????????????????????????????
    else
        indFailure=max(indFailure+1,find(nbrFailure>0,1,'first'));
    end
    [ifail,jfail] = find(MFailure(:,:,indFailure)>0);
    
    if crack_count{1,end}==0
        RVE_data.First_fragment_location=[ifail-1,jfail-1];
    end
    
    %FJ - changed to evaluate all interfaces after a crack
    % [ifail_cracks,jfail_cracks] = find(MFailure(:,:,indFailure)>0);
    
    %     % FJ - find min overlaps to determine RVE fracture toughness - this is
    %     % calculated for all strain increments
    %     MminOverlaptemp = minOverlapFracture();
    %     MminOverlap(:,:,indFailurePrevious:indFailure-1)=repmat(MminOverlaptemp(:,:,1),1,1,indFailure-indFailurePrevious);
    
    % Add cracks in material %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     AddCracks(ifail_cracks,jfail_cracks,indFailure);
    AddCracks(ifail,jfail,indFailure);
    %% FJ - recording fibre crack count
    % FJ - store cumulative crack count, indFailure and location of fibre
    % fragmentations in crack_count
    crack_count={crack_count{1,:} crack_count{1,end}+length(ifail);crack_count{2,:} indFailure;crack_count{3,:} {[(ifail-1) (jfail-1)]}};
    
    % FJ - store cumulative crack count, indFailure and location of fibre
    % fragmentations in crack_count separately for carbon and glass
    
    % reset counters
    carbon_breaks=0;
    glass_breaks=0;
    ifail_carbon=[];
    jfail_carbon=[];
    ifail_glass=[];
    jfail_glass=[];
    
    %count cracks and record crack locations  for carbon / glass only if a
    % break occurs for that fibre type
    for crack_counter=1:length(ifail)
        if RVE_outer_layer(ifail(crack_counter),jfail(crack_counter))==1
            carbon_breaks=carbon_breaks+1;
            ifail_carbon=[ifail_carbon; ifail(crack_counter)];
            jfail_carbon=[jfail_carbon; jfail(crack_counter)];
        elseif RVE_outer_layer(ifail(crack_counter),jfail(crack_counter))==0
            glass_breaks=glass_breaks+1;
            ifail_glass=[ifail_glass; ifail(crack_counter)];
            jfail_glass=[jfail_glass; jfail(crack_counter)];
        end
    end
    
    % FJ - store cumulative crack count, indFailure and location of fibre
    % fragmentations in crack_count separately for carbon and glass - comment this out for speed if you wish
        RVE_data.Carbon_crack_count={RVE_data.Carbon_crack_count{1,:} RVE_data.Carbon_crack_count{1,end}+carbon_breaks;RVE_data.Carbon_crack_count{2,:} indFailure;RVE_data.Carbon_crack_count{3,:} {[(ifail_carbon-1) (jfail_carbon-1)]}};
        RVE_data.Glass_crack_count={RVE_data.Glass_crack_count{1,:} RVE_data.Glass_crack_count{1,end}+glass_breaks;RVE_data.Glass_crack_count{2,:} indFailure;RVE_data.Glass_crack_count{3,:} {[(ifail_glass-1) (jfail_glass-1)]}};
    
    %% storing results for this iteration
    % Store results between last crack and now %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MSScurves(:,:,indFailurePrevious:indFailure-1)=Mfib(:,:,indFailurePrevious:indFailure-1);
    Mint_record(:,:,indFailurePrevious:indFailure-1)=Mint(:,:,indFailurePrevious:indFailure-1);
    % FJ - add previous locations of failure to additional failures
    MFailureOverall(:,:,indFailure:end)= MFailureOverall(:,:,indFailure:end)+repmat(MFailure(:,:,indFailure),1,1,nEpsilon-indFailure+1);
    % Show progression during simulation (not necessary) %%%%%%%%%%%%%%%%%%
    
    % FJ -
    %     C = [EpsilonGrid(indFailure) sum(sum(MSScurves(2:end-1,2:end-1,indFailure-1).*pi.*MD(2:end-1,2:end-1).^2/4))./sum(sum(pi.*MD(2:end-1,2:end-1).^2/4)).*Vf];
    
    % Break after n% to make the model quicker
    % FJ - MAKE SURE THIS IS BIGGER THAN THE MOST DUCTILE STRAIN TO
    % FAILURE!!!
    if EpsilonGrid(indFailure)>max_strain %%%%%%%%% || sum(sum(MSScurves(:,:,indFailure-1)))==0
        break;
    end
    
    1;
    
    % keep a record of the number of neighbours after each iteration -
    % comment out for speed
    %     RVE_data.MNumNeighbours_record(:,:,indFailure)=MNumNeighbour(2:end-1,2:end-1);
    
    
    
end
switch failure_modes_switch
    case 1
        Identify_failure_mode(Mint_record);
    otherwise
end

RVE_data.Cumul_crack_count=crack_count;
RVE_data.Fibre_crack_count=RVE_data.Fibre_crack_count(2:end-1,2:end-1);
RVE_data.Fibre_precrack_count=RVE_data.Fibre_precrack_count(2:end-1,2:end-1);
RVE_data.Fibre_partialfailed_count=RVE_data.Fibre_partialfailed_count(2:end-1,2:end-1);
RVE_data.Fibre_prefailed_count=RVE_data.Fibre_prefailed_count(2:end-1,2:end-1);
%
% % FJ - collect overlap lengths for each fibre
% Overlaps_top(1:n_rows,1:n_columns,:)=Moverlaps(2:2:end-2,2:1:end,:);
% Overlaps_bottom(1:n_rows,1:n_columns,:)=Moverlaps(4:2:end,2:1:end,:);
% Overlaps_left(1:n_rows,1:n_columns,:)=Moverlaps(3:2:end-1,1:1:end-1,:);
% Overlaps_right(1:n_rows,1:n_columns,:)=Moverlaps(3:2:end-1,2:1:end,:);
% % FJ - concatenate array of fibre lengths and sort smallest overlap to biggest,
% % for each fibre
% Overlaps=sort(cat(3,Overlaps_top,Overlaps_bottom,Overlaps_left,Overlaps_right),3);
% % FJ - save min and max overlap lengths for each fibre
% RVE_data.Min_overlap_lengths=Overlaps(:,:,1);
% RVE_data.Max_overlap_lengths=Overlaps(:,:,end);
%
% % capture individual cracks throughout RVE - use this for plots in
% HaNa and Marco paper
RVE_data.Total_fibre_cracks=MFailureOverall(2:end-1,2:end-1,:);

% % FJ - Plot progression of cracks throughout RVE
% figure(26)
% for ii=1:2:size(nbrFailure,2)
%     imagesc(MFailureOverall(2:end-1,2:end-1,ii));
%     pause
% end
% clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     % FJ - find min overlaps to determine RVE fracture toughness - this is
% %     % calculated for final strain increment
% MminOverlaptemp = minOverlap();
% MminOverlap=MminOverlaptemp;
% MminOverlapType=MminOverlaptemp(:,:,2);

MminOverlap=[];
MminOverlapType=[];

% Post-Process
%% %%%%%% FJ - SORAIA: this may need to change when different diameters are implemented
Stress(1:size(MSScurves,3))=sum(sum(MSScurves(2:end-1,2:end-1,:).*pi.*repmat(MD(2:end-1,2:end-1),1,1,size(MSScurves,3)).^2/4))./sum(sum(pi.*repmat(MD(2:end-1,2:end-1),1,1,size(MSScurves,3)).^2/4)).*Vf;
Strain=EpsilonGrid';
% figure (3)
% plot(Strain,Stress)
% xlim([0 3])
% hold on
Ecomposite=[Strain' Stress'];
1;

% keep a record of interface stresses - don't keep if you want this to run
% quickly!
% RVE_data.Mint_record = Mint;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OTHER SIDE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Min overlap with time (for fracture toughness calculations)

    function MminOverlaptemp = minOverlapFracture()
        
        MminOverlaptemp=zeros(size(Mint,1),size(Mint,2),2);
        
        for io=1:size(Mfib,1)-1
            for jo=1:size(Mfib,2)-1
                Icracks=findIcrack(io,jo,io+1,jo);
                min01=min((Icracks(1,2:end)-Icracks(1,1:end-1))./abs(Icracks(2,2:end)-Icracks(2,1:end-1)));
                MminOverlaptemp(2*io-1,jo,1)=min01;
                MminOverlaptemp(2*io-1,jo,2)=RVE_outer_layer(io,jo)+RVE_outer_layer(io+1,jo);
                
                Icracks=findIcrack(io,jo,io,jo+1);
                min02=min((Icracks(1,2:end)-Icracks(1,1:end-1))./abs(Icracks(2,2:end)-Icracks(2,1:end-1)));
                MminOverlaptemp(2*io,jo,1)=min02;
                MminOverlaptemp(2*io,jo,2)=RVE_outer_layer(io,jo)+RVE_outer_layer(io,jo+1);
            end
        end
        
        % removing interfaces on the edges
        MminOverlaptemp([1 2 end-1 end],[2:end-1],:)=NaN;
        MminOverlaptemp([2:end-1],1,:)=NaN;
        MminOverlaptemp([1:2:end-1],end-1,:)=NaN;
        
        fibre_overlap_record=zeros(n_rows,n_columns,4);
        % left interfaces
        fibre_overlap_record(:,:,1)=MminOverlaptemp(3:2:2*n_rows+1,1:end-2,1);
        % right interfaces
        fibre_overlap_record(:,:,2)=MminOverlaptemp(3:2:2*n_rows+1,2:end-1,1);
        % top interfaces
        fibre_overlap_record(:,:,3)=MminOverlaptemp(2:2:2*n_rows,2:end-1,1);
        % bottom interfaces
        fibre_overlap_record(:,:,4)=MminOverlaptemp(4:2:2*n_rows+2,2:end-1,1);
        
        % removing interfaces on the edges
        MminOverlaptemp([1 2 end-1 end],[2:end-1],:)=0;
        MminOverlaptemp([2:end-1],1,:)=0;
        MminOverlaptemp([1:2:end-1],end-1,:)=0;
        
        
    end

%% Min overlap with time (for fracture toughness calculations)

    function OverlapData()
        
        MminOverlaptemp=zeros(size(Mint,1),size(Mint,2),2);
        
        for io=1:size(Mfib,1)-1
            for jo=1:size(Mfib,2)-1
                Icracks=findIcrack(io,jo,io+1,jo);
                min01=min((Icracks(1,2:end)-Icracks(1,1:end-1))./abs(Icracks(2,2:end)-Icracks(2,1:end-1)));
                MminOverlaptemp(2*io-1,jo,1)=min01;
                MminOverlaptemp(2*io-1,jo,2)=RVE_outer_layer(io,jo)+RVE_outer_layer(io+1,jo);
                
                Icracks=findIcrack(io,jo,io,jo+1);
                min02=min((Icracks(1,2:end)-Icracks(1,1:end-1))./abs(Icracks(2,2:end)-Icracks(2,1:end-1)));
                MminOverlaptemp(2*io,jo,1)=min02;
                MminOverlaptemp(2*io,jo,2)=RVE_outer_layer(io,jo)+RVE_outer_layer(io,jo+1);
            end
        end
        
        % removing interfaces on the edges
        MminOverlaptemp([1 2 end-1 end],[2:end-1],:)=NaN;
        MminOverlaptemp([2:end-1],1,:)=NaN;
        MminOverlaptemp([1:2:end-1],end-1,:)=NaN;
        
        fibre_overlap_record=zeros(n_rows,n_columns,4);
        RVE_data.Fibre_overlap_record=fibre_overlap_record;
        % top interfaces
        fibre_overlap_record(:,:,1)=MminOverlaptemp(2:2:2*n_rows,2:end-1,1);
        % right interfaces
        fibre_overlap_record(:,:,2)=MminOverlaptemp(3:2:2*n_rows+1,2:end-1,1);
        % bottom interfaces
        fibre_overlap_record(:,:,3)=MminOverlaptemp(4:2:2*n_rows+2,2:end-1,1);
        % left interfaces
        fibre_overlap_record(:,:,4)=MminOverlaptemp(3:2:2*n_rows+1,1:end-2,1);
               
        RVE_data.Fibre_overlap_record(:,:,1)=fibre_overlap_record(:,:,1);
        RVE_data.Fibre_overlap_record(:,:,2)=fibre_overlap_record(:,:,2);
        RVE_data.Fibre_overlap_record(:,:,3)=fibre_overlap_record(:,:,3);
        RVE_data.Fibre_overlap_record(:,:,4)=fibre_overlap_record(:,:,4);
        
    end

%% l to n (length to number in vector)
    function n = lton(l)
        
        %n=round(l*nLib/lf+1); %OLD VERSION
        [~,n]=min(abs(lLib-l));
    end

%% n to l (number in vector to length)
    function l = ntol(n)
        l=lf*(n-1)/nLib;
    end

%% find Mpeak (max and index)
    function [Mpeakmax,indMpeakmax]=findMpeak()
        %% %%%%% FJ - Soaraia: this function may need to change so that it calculates the stress field with different diameters
        
        %         Mintm1=Mint.*repmat(Mcalc11,1,1,size(Mint,3));
        %         Mintm2=Mint.*repmat(Mcalc12,1,1,size(Mint,3));
        Mintm1=bsxfun(@times,Mint,Mcalc11); % top/left fibre
        Mintm2=bsxfun(@times,Mint,Mcalc12); % bottom/right fibre
        
        %         Minta1=Mintm1./repmat(Mcalc21,1,1,size(Mint,3));
        %         Minta2=Mintm2./repmat(Mcalc22,1,1,size(Mint,3));
        Minta1=bsxfun(@rdivide,Mintm1,Mcalc21); % top/left fibre
        Minta2=bsxfun(@rdivide,Mintm2,Mcalc22); % bottom/right fibre
        
        %         Minta1=bsxfun(@rdivide,Mintm2,Mcalc21);
        %         Minta2=bsxfun(@rdivide,Mintm1,Mcalc22);
        
        
        Mpeak(n_rows_outer_layer,n_columns_outer_layer,nEpsilon,4)=0;
        
        % FJ - original by Joel
        %         Mpeak(2:(end-1),2:(end-1),:,1)=(Mintm2(2:2:(end-2),2:(end-1),:)+Minta1(3:2:(end-1),2:(end-1),:)+Minta1(4:2:end,2:(end-1),:)+Minta2(3:2:(end-1),1:(end-2),:))./4;
        %         Mpeak(2:(end-1),2:(end-1),:,2)=(Minta2(2:2:(end-2),2:(end-1),:)+Mintm1(3:2:(end-1),2:(end-1),:)+Minta1(4:2:end,2:(end-1),:)+Minta2(3:2:(end-1),1:(end-2),:))./4;
        %         Mpeak(2:(end-1),2:(end-1),:,3)=(Minta2(2:2:(end-2),2:(end-1),:)+Minta1(3:2:(end-1),2:(end-1),:)+Mintm1(4:2:end,2:(end-1),:)+Minta2(3:2:(end-1),1:(end-2),:))./4;
        %         Mpeak(2:(end-1),2:(end-1),:,4)=(Minta2(2:2:(end-2),2:(end-1),:)+Minta1(3:2:(end-1),2:(end-1),:)+Minta1(4:2:end,2:(end-1),:)+Mintm2(3:2:(end-1),1:(end-2),:))./4;
        
        % FJ - corrected version
        Mpeak(2:(end-1),2:(end-1),:,1)=bsxfun(@rdivide,(Mintm2(2:2:(end-2),2:(end-1),:)+Minta1(3:2:(end-1),2:(end-1),:)+Minta1(4:2:end,2:(end-1),:)+Minta2(3:2:(end-1),1:(end-2),:)),MNumNeighbour(2:end-1,2:end-1));
        Mpeak(2:(end-1),2:(end-1),:,2)=bsxfun(@rdivide,(Minta2(2:2:(end-2),2:(end-1),:)+Mintm1(3:2:(end-1),2:(end-1),:)+Minta1(4:2:end,2:(end-1),:)+Minta2(3:2:(end-1),1:(end-2),:)),MNumNeighbour(2:end-1,2:end-1));
        Mpeak(2:(end-1),2:(end-1),:,3)=bsxfun(@rdivide,(Minta2(2:2:(end-2),2:(end-1),:)+Minta1(3:2:(end-1),2:(end-1),:)+Mintm1(4:2:end,2:(end-1),:)+Minta2(3:2:(end-1),1:(end-2),:)),MNumNeighbour(2:end-1,2:end-1));
        Mpeak(2:(end-1),2:(end-1),:,4)=bsxfun(@rdivide,(Minta2(2:2:(end-2),2:(end-1),:)+Minta1(3:2:(end-1),2:(end-1),:)+Minta1(4:2:end,2:(end-1),:)+Mintm2(3:2:(end-1),1:(end-2),:)),MNumNeighbour(2:end-1,2:end-1));
                
        % FJ - not used
        %         Mpeak(2:(end-1),2:(end-1),:,1)=bsxfun(@rdivide,(Mintm2(2:2:(end-2),2:(end-1),:)+Minta1(3:2:(end-1),2:(end-1),:)+Minta1(4:2:end,2:(end-1),:)+Minta2(3:2:(end-1),1:(end-2),:)),repmat(MNumNeighbour(2:(end-1),2:(end-1)),1,1,nEpsilon));
        %         Mpeak(2:(end-1),2:(end-1),:,2)=bsxfun(@rdivide,(Minta2(2:2:(end-2),2:(end-1),:)+Mintm1(3:2:(end-1),2:(end-1),:)+Minta1(4:2:end,2:(end-1),:)+Minta2(3:2:(end-1),1:(end-2),:)),repmat(MNumNeighbour(2:(end-1),2:(end-1)),1,1,nEpsilon));
        %         Mpeak(2:(end-1),2:(end-1),:,3)=bsxfun(@rdivide,(Minta2(2:2:(end-2),2:(end-1),:)+Minta1(3:2:(end-1),2:(end-1),:)+Mintm1(4:2:end,2:(end-1),:)+Minta2(3:2:(end-1),1:(end-2),:)),repmat(MNumNeighbour(2:(end-1),2:(end-1)),1,1,nEpsilon));
        %         Mpeak(2:(end-1),2:(end-1),:,4)=bsxfun(@rdivide,(Minta2(2:2:(end-2),2:(end-1),:)+Minta1(3:2:(end-1),2:(end-1),:)+Minta1(4:2:end,2:(end-1),:)+Mintm2(3:2:(end-1),1:(end-2),:)),repmat(MNumNeighbour(2:(end-1),2:(end-1)),1,1,nEpsilon));
        %
        [Mpeakmax,indMpeakmax]=max(Mpeak,[],4);
    end

%% Add crack in material
    function [] = AddCracks(ifail,jfail,indFailure)
        for kf=1:length(ifail)
            % FJ - finds which interface has the highest stress
            ifa=indMpeakmax(ifail(kf),jfail(kf),indFailure);
            if failure_modes_switch==1
                RVE_data.Fragmentation_record{ifail(kf),jfail(kf)}(ifa,2)=indFailure;
            end
            % FJ - finds location for crack based on fibre ends alignment
            % for highest loaded neighbour
            Icracks=findIcrack(ifail(kf),jfail(kf),ifail(kf)+di(ifa),jfail(kf)+dj(ifa));
            c2=find(Icracks(2,:)==2);
            Icracks(c2);
            
            kf2=1;
            maxmin=0;
            BreakPosition=0;
            while kf2<=length(c2) % Find position
                kf2temp=kf2;
                % FJ - length of overlap on left side of crack
                l1=Icracks(1,c2(kf2))-Icracks(1,c2(kf2)-1);
                if kf2<length(c2)
                    while c2(kf2)+1==c2(kf2+1)
                        kf2=kf2+1;
                        if kf2>=length(c2)
                            break;
                        end
                    end
                end
                % FJ - length of crack overlap on right side of crack
                l2=Icracks(1,c2(kf2)+1)-Icracks(1,c2(kf2));
                % FJ - find min overlap length and store it's index
                [maxmintemp,ind]=min([l1 l2]);
                
                BreakPosition=BreakPosition*(maxmintemp<maxmin)+((Icracks(1,c2(kf2temp))*(ind==2)+Icracks(1,c2(kf2))*(ind==1))*(maxmintemp>maxmin)); % crack at longer overlap
                
                kf2=kf2+1;
            end
            
            % TEST RANDOM CRACK POSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             BreakPosition=RVE(ifail(kf),jfail(kf))*rand()*lf;
            %             BreakPosition=RVE(ifail(kf),jfail(kf))+rand()*lf;
            
            
            if (BreakPosition~=0 && BreakPosition~=Icracks(1,1))
                BP=Mfail{ifail(kf),jfail(kf)};
                BP=sort([BP BreakPosition]);
                Mfail{ifail(kf),jfail(kf)}=BP;
                MstocSWei(ifail(kf),jfail(kf),1:3)=MstocSWei(ifail(kf),jfail(kf),2:4);
                MstocS(ifail(kf),jfail(kf))=MstocSWei(ifail(kf),jfail(kf),1);
                
                % REMOVE DEAD INTERFACES
                %% %%%%%% FJ - SORAIA: we may need to change the matrix of fibre diameters accordingly when number of neighbours is reduced
                MNumNeighbour(ifail(kf),jfail(kf))=max(MNumNeighbour(ifail(kf),jfail(kf))-1,1);
                MNumNeighbour(ifail(kf)+di(ifa),jfail(kf)+dj(ifa))=max(MNumNeighbour(ifail(kf)+di(ifa),jfail(kf)+dj(ifa))-1,1);
                1;
                
                RVE_data.Fibre_crack_count(ifail-1,jfail-1)=RVE_data.Fibre_crack_count(ifail-1,jfail-1)+1;
            end
        end
    end
    function PreFailedInterface(pInterface)
        
        % FJ - create list of interfaces as per normal
        ifail_int=[];
        jfail_int=[];
        % FJ - modified for cracks on outside of RVE
        for ii=2:n_rows_outer_layer-1
            % for i=1:n_rows
            % ### FJ - na modified to support naxnb rectangle RVE ###
            for jj=2:n_columns_outer_layer-1
                %     for j=1:n_columns
                ifail_int=[ifail_int, ii];
                jfail_int=[jfail_int, jj];
            end
        end
        
        iv_int=[];
        jv_int=[];
        for kij_int=1:length(ifail_int)
            ifromfail_int=[2*ifail_int(kij_int)-2,2*ifail_int(kij_int)-1,2*ifail_int(kij_int),2*ifail_int(kij_int)-1];
            jfromfail_int=[jfail_int(kij_int),jfail_int(kij_int),jfail_int(kij_int),jfail_int(kij_int)-1];
            for kijff_int=1:4
                iv_int=[iv_int ifromfail_int(kijff_int)];
                jv_int=[jv_int jfromfail_int(kijff_int)];
            end
        end
        
        listInt_def=unique([iv_int' jv_int'],'rows');
        
        %cut-out interfaces that are included in the outer layer
        listInt_def=listInt_def(listInt_def(:,1)~=2,:);
        listInt_def=listInt_def(listInt_def(:,1)~=(2*n_rows_outer_layer-2),:);
        listInt_def=listInt_def(listInt_def(:,2)~=1,:);
        %         listInt_def=listInt_def(listInt_def(:,2)~=(n_columns_outer_layer-1),:);
        
        row_remover=[];
        for ii=1:size(listInt_def,1)
            if mod(listInt_def(ii,1),2)~=0 && listInt_def(ii,2)==n_columns_outer_layer-1
                row_remover=[row_remover;ii];
            end
        end
        
        listInt_def(row_remover,:)=[];
        
        rl=rand(length(listInt_def),1);
        
        listInt_def=[rl listInt_def];
        
        Failed_interface_list=sortrows([listInt_def],1);
        Failed_interface_list=Failed_interface_list(1:floor(pInterface*length(listInt_def)),:);
        
        Failed_interface_indices=sub2ind([2*n_rows_outer_layer-1,n_columns_outer_layer],Failed_interface_list(:,2),Failed_interface_list(:,3));
        
        RVE_data.Prefailed_interface_count(Failed_interface_indices)=1;
        
        Failed_interface_fibre_index=((mod(Failed_interface_list(:,2),2)==0).*([Failed_interface_list(:,2)./2,Failed_interface_list(:,3)]))+((mod(Failed_interface_list(:,2),2)~=0).*([ceil(Failed_interface_list(:,2)./2),Failed_interface_list(:,3)]));
        Failed_interface_fibre_index=[Failed_interface_fibre_index;((mod(Failed_interface_list(:,2),2)==0).*([Failed_interface_list(:,2)./2+1,Failed_interface_list(:,3)]))+((mod(Failed_interface_list(:,2),2)~=0).*([ceil(Failed_interface_list(:,2)./2),Failed_interface_list(:,3)+1]))];
        
        for ii = 1:size(Failed_interface_fibre_index,1)
            MNumNeighbour(Failed_interface_fibre_index(ii,1),Failed_interface_fibre_index(ii,2))=MNumNeighbour(Failed_interface_fibre_index(ii,1),Failed_interface_fibre_index(ii,2))-1;
            RVE_data.Fibre_partialfailed_count(Failed_interface_fibre_index(ii,1),Failed_interface_fibre_index(ii,2))=RVE_data.Fibre_partialfailed_count(Failed_interface_fibre_index(ii,1),Failed_interface_fibre_index(ii,2))+1;
        end
        RVE_data.Prefailed_interface_count([1,2,end-1,end],:)=NaN;
        RVE_data.Prefailed_interface_count(1:end,[1,end])=NaN;
        RVE_data.Prefailed_interface_count(1:2:end,end-1)=NaN;
        
        %record interface defects for EoR heat maps
        RVE_data.Interface_partialfailed_record=zeros(n_rows,n_columns,4);
        RVE_data.Interface_partialfailed_record(:,:,1)=RVE_data.Prefailed_interface_count(2:2:end-2,2:end-1);
        RVE_data.Interface_partialfailed_record(:,:,2)=RVE_data.Prefailed_interface_count(3:2:end-2,2:end-1);
        RVE_data.Interface_partialfailed_record(:,:,3)=RVE_data.Prefailed_interface_count(4:2:end,2:end-1);
        RVE_data.Interface_partialfailed_record(:,:,4)=RVE_data.Prefailed_interface_count(3:2:end-2,1:end-2);
    end

%% pre-failed interface / fibre in material
    function [] = PreFailedFibres(pcFailed,pgFailed,Failure_type)
        % FJ - find location of carbon fibres in RVE
        [ic,jc] = find(RVE_outer_layer(2:end-1,2:end-1)==1);
        % FJ - randomly select location of pc percent of carbon fibres to
        % break
        rl=rand(length(ic),1);
        break_loc_c=sortrows([rl ic jc]);
        break_loc_c=break_loc_c(1:floor(pcFailed*length(ic)),:);
        % add one onto each value to make consistent with ifail and other
        % matrices
        break_loc_c=break_loc_c(:,2:3)+1;
        
        % FJ - find location of glass fibres in RVE
        [ig,jg] = find(RVE_outer_layer(2:end-1,2:end-1)==0);
        % FJ - randomly select location of pc percent of carbon fibres to
        % break
        rl=rand(length(ig),1);
        break_loc_g=sortrows([rl ig jg]);
        break_loc_g=break_loc_g(1:floor(pgFailed*length(ig)),:);
        % add one onto each value to make consistent with ifail and other
        % matrices
        break_loc_g=break_loc_g(:,2:3)+1;
        
        break_loc=[break_loc_c;break_loc_g];
        
        break_loc_i=break_loc(:,1);
        break_loc_j=break_loc(:,2);
        
        switch Failure_type
            case 'Full'
                Num_breaks=4;
                %             case 'Partial'
                %                 Num_breaks=1;
        end
        
        for kf=1:length(break_loc_i)
            for ifa=1:Num_breaks
                % FJ - finds location for crack based on fibre ends alignment
                % for highest loaded neighbour
                Icracks=findIcrack(break_loc_i(kf),break_loc_j(kf),break_loc_i(kf)+di(ifa),break_loc_j(kf)+dj(ifa));
                c2=find(Icracks(2,:)==2);
                Icracks(c2);
                
                kf2=1;
                maxmin=0;
                BreakPosition=0;
                while kf2<=length(c2) % Find position
                    kf2temp=kf2;
                    % FJ - length of overlap on left side of crack
                    l1=Icracks(1,c2(kf2))-Icracks(1,c2(kf2)-1);
                    if kf2<length(c2)
                        while c2(kf2)+1==c2(kf2+1)
                            kf2=kf2+1;
                            if kf2>=length(c2)
                                break;
                            end
                        end
                    end
                    % FJ - length of crack overlap on right side of crack
                    l2=Icracks(1,c2(kf2)+1)-Icracks(1,c2(kf2));
                    % FJ - find min overlap length and store it's index
                    [maxmintemp,ind]=min([l1 l2]);
                    
                    BreakPosition=BreakPosition*(maxmintemp<maxmin)+((Icracks(1,c2(kf2temp))*(ind==2)+Icracks(1,c2(kf2))*(ind==1))*(maxmintemp>maxmin)); % crack at longer overlap
                    
                    kf2=kf2+1;
                end
                
                % TEST RANDOM CRACK POSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             BreakPosition=RVE(break_loc_i(kf),break_loc_j(kf))*rand()*lf;
                %             BreakPosition=RVE(break_loc_i(kf),break_loc_j(kf))+rand()*lf;
                
                
                if (BreakPosition~=0 && BreakPosition~=Icracks(1,1))
                    BP=Mfail{break_loc_i(kf),break_loc_j(kf)};
                    BP=sort([BP BreakPosition]);
                    Mfail{break_loc_i(kf),break_loc_j(kf)}=BP;
                    MstocSWei(break_loc_i(kf),break_loc_j(kf),1:3)=MstocSWei(break_loc_i(kf),break_loc_j(kf),2:4);
                    MstocS(break_loc_i(kf),break_loc_j(kf))=MstocSWei(break_loc_i(kf),break_loc_j(kf),1);
                    
                    % REMOVE DEAD INTERFACES
                    MNumNeighbour(break_loc_i(kf),break_loc_j(kf))=max(MNumNeighbour(break_loc_i(kf),break_loc_j(kf))-1,1);
                    MNumNeighbour(break_loc_i(kf)+di(ifa),break_loc_j(kf)+dj(ifa))=max(MNumNeighbour(break_loc_i(kf)+di(ifa),break_loc_j(kf)+dj(ifa))-1,1);
                    1;
                    
                    RVE_data.Fibre_crack_count(break_loc_i(kf),break_loc_j(kf))=RVE_data.Fibre_crack_count(break_loc_i(kf),break_loc_j(kf))+1;
                end
                %                 if strcmp(Failure_type,'Partial') == 1
                %                     RVE_data.Fibre_partialfailed_count(break_loc_i(kf),break_loc_j(kf))=RVE_data.Fibre_partialfailed_count(break_loc_i(kf),break_loc_j(kf))+1;
                %                 end
            end
            
            if strcmp(Failure_type,'Full') == 1
                RVE_data.Fibre_prefailed_count(break_loc_i(kf),break_loc_j(kf))=RVE_data.Fibre_prefailed_count(break_loc_i(kf),break_loc_j(kf))+1;
            end
        end
    end
%% Add residual crack in material

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFIED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [] = ResidualCrackC(pc)
        % FJ - find location of carbon fibres in RVE
        [il,jl] = find(RVE_outer_layer==1);
        % FJ - randomly select location of pc percent of carbon fibres to
        % break
        rl=rand(length(il),1);
        ll=sortrows([rl il jl]);
        ll=ll(1:floor(pc*length(il)),:);
        ll=ll(:,2:3);
        
        for kl=1:size(ll,1)
            % break selected carbon fibres
            AddResidualCracksC(ll(kl,1),ll(kl,2));
        end
    end

    function [] = ResidualCrackG(pg)
        
        % FJ - find i and j location of glass fibres in RVE
        [il,jl] = find(RVE_outer_layer==0);
        % FJ - randomly select location of pg percent of glass fibres to
        % break
        rl=rand(length(il),1);
        ll=sortrows([rl il jl]);
        ll=ll(1:floor(pg*length(il)),:);
        ll=ll(:,2:3);
        
        for kl=1:size(ll,1)
            % add residual cracks to selected glass fibres
            AddResidualCracksG(ll(kl,1),ll(kl,2));
        end
    end

    function [] = AddResidualCracksC(irb,jrb)
        % FJ - add crack along half of fibre length - use only for HiPerDiF
        % length investigation
        if HalfLengthC==1
            BreakPosition=RandomLength(irb,jrb)+lf/2;
            %             BreakPosition=RandomLength(irb,jrb)+rand()*lf;
            % FJ - Add crack in random fibre position (normal setting)
        else
            BreakPosition=RandomLength(irb,jrb)+rand()*lf;
        end
        
        %             if (BreakPosition~=Icracks(1,1))
        
        % FJ - update failure matrix and stochastic strengths
        BP=Mfail{irb,jrb};
        BP=sort([BP BreakPosition]);
        Mfail{irb,jrb}=BP;
        
        % FJ - comment stochastic strength code below to keep strength the
        % same - only do this for HiPerDiF length investigation (when
        % length is  halved for particular fibre)
        if HalfLengthC~=1
            MstocSWei(irb,jrb,1:3)=MstocSWei(irb,jrb,2:4);
            MstocS(irb,jrb)=MstocSWei(irb,jrb,1);
            RVE_data.Fibre_crack_count(irb,jrb)=RVE_data.Fibre_crack_count(irb,jrb)+1;
            RVE_data.Fibre_precrack_count(irb,jrb)=RVE_data.Fibre_precrack_count(irb,jrb)+1;
            %             RVE_data.Fragmentation_record{irb,jrb}=cat(1,RVE_data.Fragmentation_record{irb,jrb},[NaN 0]);
        end
        
    end

    function [] = AddResidualCracksG(irb,jrb)
        % FJ - add crack along half of fibre length - use only for HiPerDiF
        % length investigation
        if HalfLengthG==1
            BreakPosition=RandomLength(irb,jrb)+lf/2;
            %             BreakPosition=RandomLength(irb,jrb)+rand()*lf;
            % FJ - Add crack in random fibre position (normal setting)
        else
            BreakPosition=RandomLength(irb,jrb)+rand()*lf;
        end
        
        %             if (BreakPosition~=Icracks(1,1))
        
        % FJ - update failure matrix and stochastic strengths
        BP=Mfail{irb,jrb};
        BP=sort([BP BreakPosition]);
        Mfail{irb,jrb}=BP;
        % FJ - comment stochastic strength code below to keep strength the
        % same - only do this for HiPerDiF length investigation (when
        % length is  halved for particular fibre)
        if HalfLengthG~=1
            MstocSWei(irb,jrb,1:3)=MstocSWei(irb,jrb,2:4);
            MstocS(irb,jrb)=MstocSWei(irb,jrb,1);
            RVE_data.Fibre_crack_count(irb,jrb)=RVE_data.Fibre_crack_count(irb,jrb)+1;
            RVE_data.Fibre_precrack_count(irb,jrb)=RVE_data.Fibre_precrack_count(irb,jrb)+1;
            %             RVE_data.Fragmentation_record{irb,jrb}=cat(1,RVE_data.Fragmentation_record{irb,jrb},[NaN 0]);
        end
        %             end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MODIFIED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Epsilon interpolation
    function [SigmaIntOnEpsilon] = EpsilonInterpolation(x,y)
        
        %         if min(y(2:end-1))==0||max(x>1E+03)
        %             1;
        %         end
        SigmaIntOnEpsilon = interp1([x ; EpsilonMax+10^(-6)],[y ; 0],EpsilonGrid);
        %         SigmaIntOnEpsilon = interp1([x ; x(end)+10^(-6)],[y ; 0],EpsilonGrid);
        
    end

%% Find list/position of cracks in an interface
    function Icracks=findIcrack(if1,jf1,if2,jf2)
        
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
    function [xInt,yInt,EInt,deInt,assert_flag]=SSInterface_EffectsOfRandomness(EE_int,DE_int,ET_int,DT_int,signDeltaK,tm_int,Sm_int,lk)
        
        if lk <= lf/512
            xInt=[0;100];
            yInt=[0;0];
            EInt=0;
            deInt=0;
            assert_flag=0;
        elseif tm_int <= 5E-04
            [EInt,deInt,xInt,yInt,assert_flag] = HBaMBiLinSP_robust_FJPostProcess_update_AvgIntervals(EE_int,signDeltaK*DE_int,ET_int,signDeltaK*DT_int,5E-04,lk,Sm_int,GiicmAvg,G,1000);
        else
            [EInt,deInt,xInt,yInt,assert_flag] = HBaMBiLinSP_robust_FJPostProcess_update_AvgIntervals(EE_int,signDeltaK*DE_int,ET_int,signDeltaK*DT_int,tm_int,lk,Sm_int,GiicmAvg,G,1000);
            %         xInt=xInt*100;
            %         [xInt,yInt,EInt,deInt] = SS_HBaM(EE_int,signDeltaK*DE_int,ET_int,signDeltaK*DT_int,tm_int,lk,Sm_int,GiicmAvg,G);
        end
        
    end
%% Combination in series
    function [xSeries,ySeries]=CombineInSeries(SSn) %{xInt,yInt,EInt,deInt,lk} in each
        
        n=length(SSn);
        
        Smaxn=[];
        
        for ki=1:n
            Smaxn=[Smaxn max(SSn{ki}{2})];
        end
        
        [Smax,iMin]=min(Smaxn);
        % changed with SP updates to effects of randomness
        %         nSmax=Smax/ds;
        nSmax=floor(Smax/ds);
        
        if nSmax==0
            xSeries=[0 EpsilonMax]';
            ySeries=[0 0]'; % big assumption here ...
        else
            
            % Loading
            u=zeros(nSmax+1,1); % Total displacement
            for kl=1:n
                u=u+SSn{kl}{1}(1:nSmax+1)/100*SSn{kl}{5}; % in mm
            end
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

%% Combination in parallel


%% plot stress with time
% FJ - does what it says
    function plotstress()
        figure
        colorbar
        % STRESS
        colormap(jet(64))
        caxis([0,5000])
        % FIBRE FAILURE
        % joelmap = [1.0000    0.2745    0.2941;0.1529    0.5725    0.9569];
        % colormap(joelmap)
        % caxis([0,1])
        axis([0.5 n_rows_outer_layer+0.5 0.5 n_columns_outer_layer+0.5])
        hold on
        for in=1:nEpsilon/3
            imagesc(MSScurves(:,:,3*in));
            pause
        end
    end

%% plot fibre distribution and overlap lengths distribution
    function plotRVE()
        %         figure(1);
        %         imagesc(RVE);
        %         figure(2);
        %         imagesc(RandomLength);
    end

%% Defines matrix of fibre-types moduli E
    function Ef = defineEf(Typef,EgAvg,EgCoV,EcAvg,EcCoV)
        
        %Typef = matrix of fibre-types. Typef=0 for G, Typef=1 for C.
        %Requires Expected (Avg) and CoV=Std/Avg values of normal distribution
        
        %randn = from N(0,1). N(Avg,Std)=Avg+N(0,1)*Std=Avg*(1+N(0,1)*CoV)
        
        Ef=(1-Typef)*EgAvg.*(1+randn(size(Typef)).*EgCoV) +...  %if G
            Typef*EcAvg.*(1+randn(size(Typef)).*EcCoV);         %if C
        
    end

%% Defines matrix of shear lag thicknesses tI
    function [tv,th] = definetI(n_rows_outer_layer,n_columns_outer_layer,Vf,nfc,dg,dc)
        
        NGridtI=1000;
        
        %Calculating average interaction thickness:
        DfAvg=(1-nfc)*dg+nfc*dc;
        Df2Avg=(1-nfc)*dg^2+nfc*dc^2;
        
        tIAvg = 1024/(405*pi) * DfAvg *...
            ( sqrt(1+ (405*pi)/1024 * Df2Avg/DfAvg^2 * (1-Vf)/Vf) - 1);
        
        %Parameter for distribution (areal density of events)
        p=1024/(81*pi^2)/tIAvg^2;
        
        %Generate random numbers
        Ftv=rand(n_rows_outer_layer-1,n_columns_outer_layer);
        Fth=rand(n_rows_outer_layer,n_columns_outer_layer-1);
        
        %Convert random numbers into random values of tI
        %Finds maximum tI required.
        %Initial searching point is 10*tIAvg, and assigns the maximum
        %for a 1.1*larger thickness than the supposed maximum
        %(to avoid extrapolation)
        tIMax=fzero(@(tI)...
            1 - (1-pi*p*tI^2/6)*erfc(sqrt(pi*p)*tI/2) ...
            - tI*sqrt(p)*exp(-pi*p*tI^2/4) ...
            -max([Ftv(:);Fth(:)]),10*tIAvg)*1.1;
        
        %Creating grid with CDF:
        %Creates grid of tI values
        tIGrid=linspace(0,tIMax,NGridtI);
        %Calculates CDF at the grid values
        CDFtIGrid = 1 - (1-pi*p*tIGrid.^2/6).*erfc(sqrt(pi*p)*tIGrid/2) ...
            - sqrt(p)*tIGrid.*exp(-pi*p*tIGrid.^2/4);
        
        %Calculating tI matrices
        tv=interp1(CDFtIGrid,tIGrid,Ftv);
        th=interp1(CDFtIGrid,tIGrid,Fth);
        
    end

%% Defines matrix of shear lag strength SI
    function [Sv,Sh] = defineSI(n_rows_outer_layer,n_columns_outer_layer,SmAvg,SmCoV)
        
        %SmV = matrix of shear strength in vertical (along i) direction.
        %SmH = matrix of shear strength in horizontal (along j) direction.
        %Requires Expected (Avg) and CoV=Std/Avg values of Weibull distribution
        
        %Weibull parameters
        SmWm=fzero(@(m) sqrt(gamma(1+2/m)/(gamma(1+1/m)^2)-1)-SmCoV, 1.2/SmCoV);
        SmWs=SmAvg/(gamma(1+1/SmWm));
        
        Sv=wblrnd(SmWs,SmWm,n_rows_outer_layer-1,n_columns_outer_layer);
        Sh=wblrnd(SmWs,SmWm,n_rows_outer_layer,n_columns_outer_layer-1);
        
    end

%% Collates interaction properties and stores them as a average or minimum value per fibre
    function Record = Record_interaction_properties(V,H,Instruction)
        
        Record=zeros(n_rows,n_rows);
        Collate=NaN(n_rows,n_columns,4);
   
        %top interaction
        Collate(2:end,:,1)=V(2:end-1,2:end-1);
        % right interaction
        Collate(:,1:end-1,2)=H(2:end-1,2:end-1);
        % bottom interaction
        Collate(1:end-1,:,3)=V(2:end-1,2:end-1);
        % left interaction
        Collate(:,2:end,4)=H(2:end-1,2:end-1);
        
        switch Instruction
            case 'average'
                Record=mean(Collate,3,'omitnan');
            case 'minimum'
                Record=min(Collate,[],3,'omitnan');
            case 'maximum'
                Record=max(Collate,[],3,'omitnan');
            case 'record'
                Record=Collate;
        end
    end

%% FJ - identiy failure mode of RVE

    function Identify_failure_mode(Mint_record)
        Failure_mode_record=cell(n_rows,n_columns);
        row_modifier=[0,1,2,1];
        column_modifier=[1,1,1,0];
        
        % for every fibre interface find the peak stress and record the
        % strain at which the peak stress was achieved
        for ii=1:n_rows
            for jj=1:n_columns
                Failure_mode_record{ii,jj}{1}=double.empty(0,2);
                for n_interface=1:4
                    peaksloc=[];

                    [~,peaksloc]=findpeaks(reshape(Mint_record(2*ii+row_modifier(n_interface),jj+column_modifier(n_interface),:),size(Mint_record,3),1));

                    if isempty(peaksloc)~=1
                        Failure_mode_record{ii,jj}{1}=cat(1,Failure_mode_record{ii,jj}{1},[peaksloc+1,reshape(Mint_record(2*ii+row_modifier(n_interface),jj+column_modifier(n_interface),peaksloc(:)+1),size(peaksloc,1),1)>5]);
                    end
                end
                Failure_mode_record{ii,jj}{1}=sortrows(Failure_mode_record{ii,jj}{1});
                [~,indPeaks]=unique(Failure_mode_record{ii,jj}{1}(:,1));
                Failure_mode_record{ii,jj}{1}=Failure_mode_record{ii,jj}{1}(indPeaks,:);
            end
        end
        % fibre failure = 1
        % fragmentation = 2
        % debonding = 3
        % softening = 4
        % neighbour fragmnetation = 5
        for ii=1:n_rows
            for jj=1:n_columns
                for kk=1:size(Failure_mode_record{ii,jj}{1},1)
        
                        if Failure_mode_record{ii,jj}{1}(kk,2)==0
                            if max(Failure_mode_record{ii,jj}{1}(kk,1)==RVE_data.Fragmentation_record{ii+1,jj+1}(:,2))==1
                                Failure_mode_record{ii,jj}{2}{kk}=2;
                            elseif Failure_mode_record{ii,jj}{1}(kk,1)==RVE_data.Fragmentation_record{ii,jj+1}(3,2) || Failure_mode_record{ii,jj}{1}(kk,1)==RVE_data.Fragmentation_record{ii+1,jj+2}(4,2) || Failure_mode_record{ii,jj}{1}(kk,1)==RVE_data.Fragmentation_record{ii+2,jj+1}(1,2) || Failure_mode_record{ii,jj}{1}(kk,1)==RVE_data.Fragmentation_record{ii+1,jj}(2,2)
%                                 Failure_mode_record{ii,jj}{2}{kk}=2;
                                % changed to avoid double counting fibre
                                % fragmentations
                                Failure_mode_record{ii,jj}{2}{kk}=5;
                            else
                                Failure_mode_record{ii,jj}{2}{kk}=3;
                            end
                        else
                            if max(Failure_mode_record{ii,jj}{1}(kk,1)==RVE_data.Fragmentation_record{ii,jj+1}(:,2))==1 || max(Failure_mode_record{ii,jj}{1}(kk,1)==RVE_data.Fragmentation_record{ii+1,jj+2}(:,2))==1 || max(Failure_mode_record{ii,jj}{1}(kk,1)==RVE_data.Fragmentation_record{ii+2,jj+1}(:,2))==1 || max(Failure_mode_record{ii,jj}{1}(kk,1)==RVE_data.Fragmentation_record{ii+1,jj}(:,2))==1
                                Failure_mode_record{ii,jj}{2}{kk}=5;
                            else
                                Failure_mode_record{ii,jj}{2}{kk}=4;
                            end
                        end
                    
                end
                if isempty(Failure_mode_record{ii,jj}{1})~=1
                    if MSScurves(ii+1,jj+1,Failure_mode_record{ii,jj}{1}(end,1))==0
                        Failure_mode_record{ii,jj}{1}=[Failure_mode_record{ii,jj}{1};[Failure_mode_record{ii,jj}{1}(end,1),0]];
                        Failure_mode_record{ii,jj}{2}=cat(2,Failure_mode_record{ii,jj}{2},1);
                    end
                end
            end
        end
        
        RVE_data.Failure_mode_record=Failure_mode_record;
        
    end



end