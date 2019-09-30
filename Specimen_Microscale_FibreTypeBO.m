

%% FJ - HiperDiF wrapper function
% Specimen_dataset = Specimen_Microscale_FibreTypeBO(Vf,Carbon_fibre_name,Glass_fibre_name,Vc,lf,SmAvg,G,GiicmAvg,SmCoV,Tu0,lspecimen,n_rows0,n_columns0,pc,pInterface,pcFailed,pg,pgPartial,pgFailed,misalignment_switch,rng_switch,failure_modes_switch,parallelisation_flag)
% FJ - taken from "PostProcess_IntFullSS" script by Joël Henry. Adapted for
% support of HiPerDiF experiments. Front end code to call other functions
% for modelling of HiPerDiF specimens

function Specimen_dataset = Specimen_Microscale_FibreTypeBO(single_multi_flag,objective_fun_flag,Vf,Carbon_fibre_name,Glass_fibre_name,Vc,lf,SmAvg,G,GiicmAvg,SmCoV,Tu0,lspecimen,n_rows0,n_columns0,pc,pInterface,pcFailed,pg,pgPartial,pgFailed,misalignment_switch,rng_switch,failure_modes_switch,parallelisation_flag)

if max(abs([pc,pInterface,pcFailed,pg,pgFailed]))>0
    defect_case='defective';
else
    defect_case='pristine';
end

repeat_counter=1;
%(repeat_counter,Composite_index (1 for HMC_EG, 2 for HSC), lf(in microns),
%rest in percent
dbstop if error
lspecimen=lspecimen/1000;
lf=lf/1000;
pc=pc/100;
pInterface=pInterface/100;
pcFailed=pcFailed/100;
pg=pg/100;
pgPartial=pgPartial/100; %legacy
pgFailed=pgFailed/100;

libGen=2;

% try
%     error('error made');
% catch e
%     1;
% end

switch rng_switch
    case 0
        rng(RVE_loop_counter);
    case 1
        
        rng('shuffle','simdTwister')
    otherwise
        rng(rng_switch);
end

% rng(repeat_counter);
% switch Composite_index
%     case 1
%         % high modulus/E-glass hybrid with epoxy matrix
%         Composite_type='HMC_EG_BO';
%     case 2
%         % high strength carbon with epoxy matrix
%         Composite_type='HSC_BO';
%     case 3
%         % High strength carbon with polypropylene matrix
%         Composite_type='HSC_PP_BO';
%     case 4
%         % Joel's inputs for his heat maps {Henry2017}
%         Composite_type='HSC_Joel_BO';
%     case 5
%         % same as HMC/EG but with double strength glass
%         Composite_type='HMC_strongEG_BO';
% end
% [ec,EcAvg,EcCoV,eg,EgAvg,EgCoV,G,GiicmAvg,gmAvg,lcWei,lgWei,mcWei,mgWei,Sc,Sg,SmAvg,SmCoV,Tu0,Vf]=load_inputs_BO(Composite_type);

switch misalignment_switch
    case 0 % no misalignment
        Alignment_method='None'; % 'None', 'Misalignment', or 'Realignment'
        Alignment_distribution='Normal_distribution';
        Sigma_alignment=[0];
        
    case 1 % HiPerDiF misalignemnt
        Alignment_method='Misalignment'; % 'None', 'Misalignment', or 'Realignment'
        Alignment_distribution='Yu2014';
        Sigma_alignment=[0];
        
    case 2 % poor misalignment {Sanadi1985}
        Alignment_method='Misalignment'; % 'None', 'Misalignment', or 'Realignment'
        Alignment_distribution='Sanadi1985';
        Sigma_alignment=[0];
        
    case 3 % HiPerDiF misalignemnt with realignment
        Alignment_method='Realignment'; % 'None', 'Misalignment', or 'Realignment'
        Alignment_distribution='Yu2014';
        Sigma_alignment=[0];
        
    case 4 % poor misalignment {Sanadi1985} with realignment
        Alignment_method='Realignment'; % 'None', 'Misalignment', or 'Realignment'
        Alignment_distribution='Sanadi1985';
        Sigma_alignment=[0];
        
    case 5 % poor misalignment - Normal distro
        Alignment_method='Misalignment'; % 'None', 'Misalignment', or 'Realignment'
        Alignment_distribution='Normal_distribution';
        Sigma_alignment=[20];
        
    case 6 % good misalignment - Normal distro
        Alignment_method='Misalignment'; % 'None', 'Misalignment', or 'Realignment'
        Alignment_distribution='Normal_distribution';
        Sigma_alignment=[5];
        
    otherwise 
        Alignment_method='Misalignment'; % 'None', 'Misalignment', or 'Realignment'
        Alignment_distribution='Normal_distribution';
        Sigma_alignment=[misalignment_switch];
end

dbstop if error
if strcmp('MO_variability',single_multi_flag)==1 || strcmp('lf_test',single_multi_flag)==1 || strcmp('general',single_multi_flag)==1
    filename='.mat';
%     objective_fun_flag=char(objective_fun_flag);
else
    tmpname=tempname;
    slash_occurances=strfind(tmpname,filesep);
    filename=tmpname(slash_occurances(end)+1:end);
end

deltaEpsilon=0.02;
deltaSigma=3;

Carbon_fibre=Fibre(Carbon_fibre_name);
ec=Carbon_fibre.ef;
EcAvg=Carbon_fibre.EfAvg;
EcCoV=Carbon_fibre.EfCoV;
lcWei=Carbon_fibre.lfWei;
mcWei=Carbon_fibre.mfWei;
Sc=Carbon_fibre.Sf;
dc=Carbon_fibre.df;

Glass_fibre=Fibre(Glass_fibre_name);
eg=Glass_fibre.ef;
EgAvg=Glass_fibre.EfAvg;
EgCoV=Glass_fibre.EfCoV;
lgWei=Glass_fibre.lfWei;
mgWei=Glass_fibre.mfWei;
Sg=Glass_fibre.Sf;
dg=Glass_fibre.df;

gmAvg=SmAvg/G;

% %% Material properties
% G=1500;
%
% % frictional stress component of fibre pull-out {Zhonga}
% Tu0=[10];
% GiicmAvg=0.8;
% SmAvg=60;
% SmCoV=0.15;
% %FJ - fracture toughness of glass-glass interface
% % Giicm_GG=[GiicmAvg];
% gmAvg=SmAvg/G;
% % % FJ - Strength of glass-glass matrix interface
% % SmAvg=[SmAvg];
%
% % % FJ - max strain of glass-glass matrix interface
% % gm_GG=SmAvg/G;
%
% % %FJ - fracture toughness of glass-carbon interface
% % Giicm_GC=[GiicmAvg];
% % % FJ - Strength of glass-carbon matrix interface
% SmAvg=[80];
% Tu0=4;
% gmAvg=SmAvg/G;
% Vf=30/36*0.30;
% % % FJ - max strain of glass-carbon matrix interface
% % gm_GC=SmAvg/G;
%
% % %FJ - fracture toughness of carbon-carbon interface
% % Giicm_CC=[GiicmAvg];
% % % FJ - Strength of carbon-carbon matrix interface
% % SmAvg=[SmAvg];
% %
% % % FJ - max strain of carbon-carbon matrix interface
% % gm_CC=[SmAvg/G];
%
% % FJ - fibre volume fraction (glass and carbon combined)
% Vf=[0.30];
%
% Vc=[49/149];
%
% mcWei=[5];
%
% mgWei=[5.7119];
%
% lcWei=1;
%
% lgWei=25;
%
% Sc=3430;
%
% % FJ - diameter of carbon fibres [mm]
% dc=0.01;
% % dc=0.007;
% % FJ - stiffness of carbon fibres [N/mm^2]
% EcAvg=860000;
% EcCoV=0.05;
% % Ec=225000;
% % FJ - strain to failure of carbon fibres [%]
% ec=Sc/EcAvg*100;
% Sg=[2400];
% % Sg=[4344];
% % FJ - diameter of glass fibres [mm]
% dg=0.007;
% % dg=0.007;
% % FJ - stiffness of glass fibres [N/mm^2]
% EgAvg=73000;
% EgCoV=0.05;
% % Eg=225000;
% % FJ - strain to failure of glass fibres [%]
% eg=Sg/EgAvg*100;



%% geometric properties

% HiPerDiF is 18x387
% n_rows0=[80];
% n_columns0=[80];
c_width=[6];
n_columns_c=[2];

% FJ - define RVE configuration using x_section:
% x_section_arrangement=1: checkerboard
% x_section_arrangement=2: intermingled sandwich
% x_section_arrangement=3: intermingled
% x_section_arrangement=4: single material
% x_section_arrangement=5: HiPerDiF arrangement - not used
% x_section_arrangement=7: Dispersed 3
% x_section_arrangement=8: DI2+
% x_section_arrangement=9: blocked
% x_section_arrangement=10: DI1
% x_section_arrangement=13: checkerboard with squares size c_width
% x_section_arrangement=14: "Prefectly spaced"
% x_section_arrangement=16: "RVE_spreadsheet.mat"
% x_section_arrangement=17: 4-blocks
% x_section_arrangement=18: 9-blocks
% x_section_arrangement=25: 16-blocks
% x_section_arrangement=26: 25-blocks
% x_section_arrangement=19: Thin ply (single layer of carbons)
% x_section_arrangement=27: all carbon fibres
% FJ - define as integer
% x_section arrangement for intermingled, DI1, DI2, DI3, blocked
x_section_arrangement=int8([3]);

% FJ - set accuracy of fibre arrangement (fuzzyness)
% HiPerDiF_fuzzyness=0: nofuzzyness
% HiPerDiF_fuzzyness=1: minimal fuzzyness to enable correct nc
% HiPerDiF_fuzzyness=2: "fuzzyness fraction" - randomized
% HiPerDiF_fuzzyness=3: "fuzzyness fraction" - more controlled - includes
% correct nc
% HiPerDiF_fuzzyness=4: row-specific "fuzzy fraction"
% HiPerDiF_fuzzyness=5: block-specific "fuzzy fraction"
% HiPerDiF_fuzzyness=6: "fuzzy fraction" with small random movement from
% HiPerDiF_fuzzyness=9: "fuzzy fraction" with normal distribution of fibres
% and random offset of normal distribution mean through along width of
% sample - gives better dispersion at extremities and strata-like effect
% original row
HiPerDiF_fuzziness=int8([1]);


% define amount of fuzziness in cross-section arrangement
level_of_fuzziness=[0.8];

% FJ - set fibre overlap type
% overlap_type==1: standard random overlap
% overlap_type==2: HiPerDiF overlap + tolerance
% overlap_type==3: perfectly staggered
overlap_type=int8(1);
% % RVE length
% lf=[3];
% % Specimen length
% lspecimen=lf;

% FJ - calculate the number of representative volume elements along
% specimen length
nRVE=floor(lspecimen/lf);

%% sum of modulus, thickness, etc.
% FJ - calculate equivalent platelet thickness on one side of fibre
% boundary
Tc=dc/4; % equivalent thickness of carbon
Tg=dg/4; % equivalent thickness of glass

% calculate matrix thickness for library generation
t1=[0 0;1 1];
t2=[0 1;0 1];
% FJ - calculate sigma(Ea, Eb) for different  fibre types (c-c,c-g;c-g,g-g)
EEm=EcAvg*(t1+t2)+EgAvg*(2-t1-t2);
% FJ - calculate delta(Ea, Eb) for different  fibre types (c-c,c-g;c-g,g-g)
DEm=(EcAvg-EgAvg)*(t1-t2);
% FJ - calculate sigma(Ta, Tb) for different  fibre types (c-c,c-g;c-g,g-g)
ETm=Tc*(t1+t2)+Tg*(2-t1-t2);
% FJ - calculate delta(Ta, Tb) for different  fibre types (c-c,c-g;c-g,g-g)
DTm=(Tc-Tg)*(t1-t2);

% FJ - calculate glass fibre volume fraction (relative to # carbon fibres)
Vg=1-Vc;
% ### FJ - na modified to support naxnb rectangle RVE ###
nf=n_rows0.*n_columns0;        % Total number of fibres

% FJ - calculate number of carbon fibres, based on area fraction
nc=floor((n_rows0*n_columns0)*(Vc/dc^2)/(Vc/dc^2+(1-Vc)/dg^2));

% FJ - calculate number of glass fibres
ng=nf-nc;

if pcFailed>0 || pgFailed>0
        nc0=nc;
        ng0=ng;
        nf0=nf;
        Vf0=Vf;
        Vc0=Vc;
        AR0=n_columns0/n_rows0; % check aspect ratio %Eqn 13 {Finley2019}


        %             Vd=Vf0*(pcFailed*nc0*dc^2+pgFailed*ng0*dg^2)/(nc0*dc^2+ng0*dg^2);

        %             Vf2=Vf0*(Vf0/(Vf0-Vd));

        %             Vc2=Vc0^2*(1+((ng0*dg^2*(1-pgFailed))/(nc0*dc^2*(1-pcFailed))));

        nc=nc0*(1+pcFailed);
        ng=ng0*(1+pgFailed);
        nf=nc+ng;

        n_rows=round(sqrt((nc0*(1+pcFailed)+ng0*(1+pgFailed))*(n_rows0/n_columns0))); %Eqns 14 and 15 {Finley2019}
        n_columns=round(sqrt((nc0*(1+pcFailed)+ng0*(1+pgFailed))*(n_columns0/n_rows0)));

        AR=n_columns/n_rows; %Eqn 13 {Finley2019}

        Vf=Vf*((nc0*(1+pcFailed)*dc^2+ng0*(1+pgFailed)*dg^2)/(nc0*(1)*dc^2+ng0*(1)*dg^2)); %Eqn 16 {Finley2019}
        Vc=(nc0*(1+pcFailed)*dc^2)/(nc0*(1+pcFailed)*dc^2+ng0*(1+pgFailed)*dg^2); %Eqn 17 {Finley2019}
else
    n_rows=n_rows0;
    n_columns=n_columns0;
end

% calculate fuzzy factor, based on nc and level of fuzzyness
Fuzzy_Factor=((2*nc)/(n_rows*n_columns)*(1-Vc)).*level_of_fuzziness;

% Generation of the matrix thickness (used for library generation only)
% FJ - equation (2) in ref [1], although it appears to be diffferent?
tm=(sqrt(-Vf*(-((nc/nf)+(ng/nf))*((nc/nf)*dc^2+(ng/nf)*dg^2)*pi+(nc/nf)*Vf*(ng/nf)*(dc-dg)^2*pi)*pi)+(-(nc/nf)*dc-(ng/nf)*dg)*Vf*pi)/(Vf*pi*((nc/nf)+(ng/nf)));

%% Debugging parameters - not needed unless somehting is fishy
% lf_SPSL_const = [1/8 1/4 1/2 1];
lf_SPSL = [1];
% lf_hybrid_const=[1/8 1/4 1/2 1];
lf_hybrid=[1];
Sc_frag_const=[1];

RVEg=1;

%% Library generation parameters

% library generation switch
%LigGen=0: load library
%LigGen=1: create library for use on these specimens
%LigGen=2: evaluate all interfaces individually

switch libGen
    case 0
        load
    case 1
        % %FJ - specify the number of different overlap lengths the library will
        % %generate stress-strain curves for
        nLib=128; % number of different lengths available
        %  % FJ - call libraryGeneration function to generate libraary data
        [AA,AAE,AAde, AB,ABE,ABde,BB,BBE,BBde,deltaSigma,lLib]=libraryGeneration(EEm,DEm,ETm,DTm,tm,lf,nLib,GiicmAvg,SmAvg,gmAvg,GiicmAvg,SmAvg,gmAvg,GiicmAvg,SmAvg,gmAvg,lf_SPSL,lf_hybrid);
    case 2
        AA=[];AAE=[]; AAde=[]; AB=[]; ABE=[];ABde=[];BB=[];BBE=[];BBde=[];lLib=[];
end

%% residual crack parameters - unused for effects of randomness study
% add p (%) of crack in fibres
% pc=0.0;
% pg=0.0;

% HiPerDiF length study
HalfLengthC=0;
HalfLengthG=0;
% HalfLengthC=[0 0 1 0];
% HalfLengthG=[0 1 0 0];




%% stress interpolation settings
% interpolation method - choose how to interpolate stresses to form
% specimen from RVEs 1 = hybrid interpolation (trimmed to fracture strain
% first), 2= end-to-start interplation, 3 = start-to-end interpolation, 4 =
% interpolation on strain, 5 = intepolation on strain with smart cut-off.
interpolation_method=int8([5]);

for i=1:length(interpolation_method)
    if interpolation_method(i)~=1 && interpolation_method(i)~=2 && interpolation_method(i)~=3 && interpolation_method(i)~=4 && interpolation_method(i)~=5
        error('incorrect interpolation method selection');
    end
end

%% Fracture method
% Fracture_method =1: Joel's, =2: Fin's
Fracture_method = [2];

%% main code for specimen-level analysis

if overlap_type==1
    % FJ - fibres are same length but are staggered at random intervals - set
    % fibre length at intersection with start of RVE
    % keep random seed same for each run
    %     rng(1);
    % random overlap
    first_fibre_overlap=rand((n_rows+2),(n_columns+2))*lf; % for  same end-fibre locations
    
elseif overlap_type==2
    % unused for effects of randomnes study
    overlap_range=0.2;
    % Fj - apply overlap tolerance to give first fibre length +/- overlap
    % tolerance - means that first fibres may be a bit longer than defined?
    first_fibre_overlap=rand((n_rows+2),(n_columns+2))*overlap_range+(lf-overlap_range/2);
    
elseif overlap_type==3
    first_fibre_overlap=zeros((n_rows+2),(n_columns+2));
    
    first_fibre_overlap(1:2:(n_rows+2),1:2:(n_columns+2))=1;
    first_fibre_overlap(2:2:(n_rows+2),2:2:(n_columns+2))=1;
    first_fibre_overlap=first_fibre_overlap*lf/2;
end


% set max strain as 2* max strain to failure of fibres (and douible it if
% Weibull refernce strength is large (suggests a stronger material)
    max_strain=max(2*max(ec)+((lcWei>25)*(lcWei/25-1)*2*max(ec)),2*max(eg)+((lgWei>25)*(lgWei/25-1)*2*max(eg)));


% RVE_data(1:nRVE)=RVE_dataset;

Specimen_dataset=Specimen;

% FJ - loop through each RVE calcaulating the stress-strain response
if parallelisation_flag==1 && floor(lspecimen/lf)>1
    parfor i=1:nRVE
        % FJ - call the intermingled model to generate the stress-strain
        % response for each RVE and store all results in a cell array
        %     [EE{i},MSScurves{i},alphag(i),alphac(i),MminOverlap{i},MminOverlapType{i}]=HiPerDiF_ModelFullSSStrainStep_RectangleXSection(Vf,lf,pc,pg,HalfLengthC,HalfLengthG,dg,Eg,Sg,dc,Ec,Sc,Vc,G,RVEg,first_fibre_overlap,n_rows,n_columns,HiPerDiF_fuzzyness,Fuzzy_Factor,x_section_arrangement,max_strain,Giicm_GG,SmAvg,gm_GG,Giicm_GC,SmAvg,gm_GC,Giicm_CC,SmAvg,gm_CC,c_width,n_columns_c,mcWei,mgWei,lcWei,lgWei,Sc_frag_const,lf_SPSL,lf_hybrid,...
        %         AA,AAE,AAde, AB,ABE,ABde, BB,BBE,BBde, deltaSigma, lLib);
        [EE{i},MSScurves{i},alphag(i),alphac(i),MminOverlap{i},MminOverlapType{i},RVE_data{i},tm_temp(i)]=RVE_Microscale_BO(i*repeat_counter,deltaEpsilon,Vf,lf,pc,pInterface,pcFailed,pg,pgPartial,pgFailed,HalfLengthC,HalfLengthG,dg,EgAvg,EgCoV,Sg,dc,EcAvg,EcCoV,Sc,Vc,G,RVEg,first_fibre_overlap,n_rows,n_columns,HiPerDiF_fuzziness,Fuzzy_Factor,x_section_arrangement,max_strain,GiicmAvg,SmAvg,SmCoV,gmAvg,c_width,n_columns_c,mcWei,mgWei,lcWei,lgWei,Sc_frag_const,lf_SPSL,lf_hybrid,libGen,AA,AAE,AAde, AB,ABE,ABde, BB,BBE,BBde, deltaSigma, lLib,rng_switch,failure_modes_switch);
        %EE{i}=HiPerDiF_ModelFullSSStrainStep(Vf,lf,dg,Eg,Sg,dc,Ec,Ec,Vc,G,RVEg,first_fibre_overlap,x_section_arrangement)
        % FJ - loop counter
        i
    end
else
    for i=1:nRVE
        % FJ - call the intermingled model to generate the stress-strain
        % response for each RVE and store all results in a cell array
        %     [EE{i},MSScurves{i},alphag(i),alphac(i),MminOverlap{i},MminOverlapType{i}]=HiPerDiF_ModelFullSSStrainStep_RectangleXSection(Vf,lf,pc,pg,HalfLengthC,HalfLengthG,dg,Eg,Sg,dc,Ec,Sc,Vc,G,RVEg,first_fibre_overlap,n_rows,n_columns,HiPerDiF_fuzzyness,Fuzzy_Factor,x_section_arrangement,max_strain,Giicm_GG,SmAvg,gm_GG,Giicm_GC,SmAvg,gm_GC,Giicm_CC,SmAvg,gm_CC,c_width,n_columns_c,mcWei,mgWei,lcWei,lgWei,Sc_frag_const,lf_SPSL,lf_hybrid,...
        %         AA,AAE,AAde, AB,ABE,ABde, BB,BBE,BBde, deltaSigma, lLib);
        [EE{i},MSScurves{i},alphag(i),alphac(i),MminOverlap{i},MminOverlapType{i},RVE_data{i},tm_temp(i)]=RVE_Microscale_BO(i*repeat_counter,deltaEpsilon,Vf,lf,pc,pInterface,pcFailed,pg,pgPartial,pgFailed,HalfLengthC,HalfLengthG,dg,EgAvg,EgCoV,Sg,dc,EcAvg,EcCoV,Sc,Vc,G,RVEg,first_fibre_overlap,n_rows,n_columns,HiPerDiF_fuzziness,Fuzzy_Factor,x_section_arrangement,max_strain,GiicmAvg,SmAvg,SmCoV,gmAvg,c_width,n_columns_c,mcWei,mgWei,lcWei,lgWei,Sc_frag_const,lf_SPSL,lf_hybrid,libGen,AA,AAE,AAde, AB,ABE,ABde, BB,BBE,BBde, deltaSigma, lLib,rng_switch,failure_modes_switch);
        %EE{i}=HiPerDiF_ModelFullSSStrainStep(Vf,lf,dg,Eg,Sg,dc,Ec,Ec,Vc,G,RVEg,first_fibre_overlap,x_section_arrangement)
        % FJ - loop counter
        i
    end
end




tm=mean(tm_temp);
%
% save(['Experiment_data_' num2str(repeat_counter)],'EE','MSScurves','alphag','alphac','crack_count','RVE_data','-v7.3');
%
% load('Experiment_data_master_150.mat')

for i=1:nRVE
    Specimen_dataset.RVE{i}=RVE_data{i};
end

% FJ - clear RVE_data to increase available RAM
clear RVE_data;

%% FJ - calculate stress-strain response of RVE chain over complete specimen length
% FJ - find the max strength of the strongest RVE
maxi=0;
for i=1:nRVE
    maxi=max(maxi,max(EE{i}(:,2)));
end
maxi;

% FJ - find the strength of the weakest RVE
mini=maxi;
for i=1:nRVE
    mini=min(mini,max(EE{i}(:,2)));
end
mini;

for i=1:nRVE
    [maxX(i) e_maxX(i)]=max(EE{i}(:,2));
end

[min_maxX i_min_maxX]=min(maxX);

% FJ - define the stress increment for plots
deltaSigma=3;
% FJ - set strength scale using highest strength value from all RVEs
% FJ - check whether this should be based on max RVE strength or min RVE
% strength. Why would it be different between stress and strain?
NewSigmaGrid=[0:deltaSigma:maxi+deltaSigma]';
% FJ - calculate number of strain points by dividing the weakest RVE by the
% strength increment and then round down to closest integer
% neps=floor(mini/deltaSigma);
neps=e_maxX(i_min_maxX);
% FJ - pre-allocate Eps outside of loop
Eps=[];
Specimen_SS_Response_unprocessed_temp=cell(nRVE,1);

% FJ - interpolate strains to give RVE stress-strain response up to
% lowest RE strength with common stress-strain scale across RVEs
for i=1:nRVE
    % FJ - define x and y as a vector of RVE strains and RVE strengths respectively - this is overwritten with
    % each loop
    x=EE{i}(:,1);
    y=EE{i}(:,2);
    
    % FJ - find the max strength of the current RVE and the index where
    % this occurs
    [m,ind]=max(y);
    
    % FJ - output unprocessed stress-strain curve
    x_unprocessed=x;
    y_unprocessed=y;
    
    Specimen_dataset.RVE{i}.Stress_strain_unprocessed=[x(1:neps),y(1:neps)];
    
    
    % FJ - trim the contents of the matrices x and y up until the strain at
    % which the weakst RVE is at its max
    x=x(1:neps);
    y=y(1:neps);
    
    
    
    
    
    if interpolation_method == 3
        % start-to-end interpolation FJ - Get rid of load drops in stress-strain curve - why is this done? to help with interpolation?
        % Larger epsilon
        n=1;
        while n<length(x)-1
            if y(n+1)<y(n)
                y=[y(1:n); y(n+2:end)];
                x=[x(1:n); x(n+2:end)];
            else
                n=n+1;
            end
        end
        
        % FJ - interpolate data using common scales for all RVEs
        %     NewEpsilonGrid = x;
        NewEpsilonGrid = interp1(y,x,NewSigmaGrid);
        %         NewEpsilonGrid = interpn(y,x,NewSigmaGrid,'cubic');
        % FJ - trim the strain vector to the number of strain points
        NewEpsilonGrid=NewEpsilonGrid(1:neps);
        
        
        
        % FJ - populate matrix of strains
        Eps(:,i)=NewEpsilonGrid;
        
    elseif interpolation_method == 1 || interpolation_method == 2
        %Smaller epsilon
        n=0;
        while n<length(x)-1
            if y(end-(n+1))>y(end-n)
                y=[y(1:end-(n+1)-1); y(end-(n+1)+1:end)];
                x=[x(1:end-(n+1)-1); x(end-(n+1)+1:end)];
            else
                n=n+1;
            end
        end
        
        % FJ - interpolate data using common scales for all RVEs
        %     NewEpsilonGrid = x;
        NewEpsilonGrid = interp1(y,x,NewSigmaGrid);
        disp([NewEpsilonGrid,NewSigmaGrid]);
        %         NewEpsilonGrid = interpn(y,x,NewSigmaGrid,'cubic');
        % FJ - trim the strain vector to the number of strain points
        NewEpsilonGrid=NewEpsilonGrid(1:neps);
        
        
        
        % FJ - populate matrix of strains
        Eps(:,i)=NewEpsilonGrid;
        
    else
        Strengths(:,i)=y(1:neps);
    end
    
    
end

if interpolation_method == 4 || interpolation_method == 5
    EpsSeries=x(1:neps);
    SigSeries=sum(Strengths,2)/nRVE;
else
    % FJ - calculate the average strain vector across all RVEs
    EpsSeries=sum(Eps,2)/nRVE;
    
    % FJ - trim the strength vector so that the max strength is that of the
    % weakest RVE (chain theory)
    SigSeries=NewSigmaGrid(1:neps);
    
end

Specimen_SS_Response1_unprocessed(:,1)=x_unprocessed(1:neps);
Specimen_SS_Response1_unprocessed(:,2)=y_unprocessed(1:neps);
% Specimen_data{:}.Stress_strain_unprocessed=Specimen_SS_Response1_unprocessed;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Estimating the fracture toughness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Enom=SigSeries(2)/EpsSeries(2);

lcWei=lcWei;   mcWei=mcWei;        XcWei=Sc; Dfc=dc;  % Carbon
lgWei=lgWei;   mgWei=mgWei;       XgWei=Sg; Dfg=dg;  % Glass


% Gm=mean([Giicm_GG Giicm_GC Giicm_CC]);
Gm=GiicmAvg;
% SmAvg=mean([SmAvg(end) SmAvg(end) SmAvg(end)]);


% Tu=10*(1-Vf);           % Matrix
% Tu0=10;

% RSc=(Vc/Dfc)/(Vc/Dfc+(1-Vc)/Dfg);
% RSg=((1-Vc)/Dfg)/(Vc/Dfc+(1-Vc)/Dfg);

RSc=Vc;
RSg=1-Vc;
% Df=Dfc*RSc+Dfg*RSg;
Df=sqrt(Dfc^2*RSc+Dfg^2*RSg);

lfG=(0.1:0.05:(lf*1.33))';
alphag=mean(alphag); %0.4518
alphac=mean(alphac); %0.4140
% alphac=1;
% alphag=1;

%% method 1 for fracture toughness - from ICCM19 (soraia)
% % Carbon
% % FJ - l_o in eqn 22 {S.Pimenta}
% lemb0c=XcWei*Dfc/(4*max(Sm))*(lcWei./lfG.*alphac.*(mcWei+1)).^(1/mcWei);
% % lemb0c=XcWei*Dfc/(4*max(Sm))*(lcWei./lfG.*alphac).^(1/mcWei);
% % FJ - eq.22 in {S.Pimenta}
% Gdebc=4*Vf*Gm./(Dfc*lfG).*lemb0c.^2/mcWei.*(gammainc((lfG/2./lemb0c).^mcWei,2/mcWei)*gamma(2/mcWei));
% Gpoc =2*Vf*Tu./(Dfc*lfG).*lemb0c.^3/mcWei.*(gammainc((lfG/2./lemb0c).^mcWei,3/mcWei)*gamma(3/mcWei));
% Gcrit_c=Gdebc+Gpoc;
% Gcrit_c(find(Gcrit_c==max(Gcrit_c)):end)=max(Gcrit_c); %remove decrease for long overlap
% Gcrit_c=interp1q(lfG,Gcrit_c,lf);
% Gcrit_c=repmat(Gcrit_c,1,nRVE);
%
% % Glass
% % FJ - l_o in eqn 22 {S.Pimenta}
% lemb0g=XgWei*Dfg/(4*max(Sm))*(lgWei./lfG.*alphag.*(mgWei+1)).^(1/mgWei);
% % lemb0g=XgWei*Dfg/(4*max(Sm))*(lgWei./lfG.*alphag).^(1/mgWei);
% % FJ - eq.22 in {S.Pimenta} - fracture toughness normalised by debonding
% % area
% Gdebg=4*Vf*Gm./(Dfg*lfG).*lemb0g.^2/mgWei.*(gammainc((lfG/2./lemb0g).^mgWei,2/mgWei)*gamma(2/mgWei));
% Gpog =2*Vf*Tu./(Dfg*lfG).*lemb0g.^3/mgWei.*(gammainc((lfG/2./lemb0g).^mgWei,3/mgWei)*gamma(3/mgWei));
% Gcrit_g=Gdebg+Gpog;
% Gcrit_g(find(Gcrit_g==max(Gcrit_g)):end)=max(Gcrit_g); %remove decrease for long overlap
% Gcrit_g=interp1q(lfG,Gcrit_g,lf);
% Gcrit_g=repmat(Gcrit_g,1,nRVE);

%
% % Composite - FJ - commented out for variable Gcrit throughout cross-section
% Gcrit=RSc*(Gdebc+Gpoc)+RSg*(Gdebg+Gpog);
% Gcrit(find(Gcrit==max(Gcrit)):end)=max(Gcrit); %remove decrease for long overlap
% Gcrit=interp1q(lfG,Gcrit,lf);
% Gcrit=repmat(Gcrit,1,nRVE);

% Gcrit_c=repmat(Gcrit_c,1,nRVE);
% Gcrit_g=repmat(Gcrit_g,1,nRVE);

%% method 2 for fracture toughness (Joel's pull-out) - remember to extend max strain until Wtotal settles!
% for i=1:nRVE
% %     pause
%
%     % FJ - retained to echeck whether minimum Wtotal is being reached
%     diam=(MminOverlapType{i}==0).*Dfg+(MminOverlapType{i}==1).*(Dfg+Dfc)/2+(MminOverlapType{i}==2).*Dfc;
%     area=(MminOverlapType{i}==0).*(pi*(Dfg/2)^2/2/Vf)+(MminOverlapType{i}==1).*((pi*(Dfg/2)^2/2+pi*(Dfc/2)^2/2)/2/Vf)+(MminOverlapType{i}==2).*(pi*(Dfc/2)^2/2/Vf);
%     Wdeb=bsxfun(@times,pi.*diam/4*Gm,MminOverlap{i});
%     Wpo=bsxfun(@times,pi.*diam/4*Tu,MminOverlap{i}.^2/2);
%     Wtotal(1:size(MminOverlap{i},3),i)=sum(sum(Wdeb+Wpo))/sum(sum(area));
%
% %     pause
%
%     % FJ - calculate diameter of carbons
%     diam_c=bsxfun(@times,MminOverlapType{i}==2,Dfc);
%     % Calculate area of carbons
%     area_c=bsxfun(@times,MminOverlapType{i}==2,(pi*(Dfc/2)^2/2/Vf));
%     % Calculate work for debonding of carbons
%     Wdeb_c=bsxfun(@times,pi.*diam_c/4*Gm,MminOverlap{i});
%     % calculate work for pull-out of carbons
%     Wpo_c=bsxfun(@times,pi.*diam_c/4*Tu,MminOverlap{i}.^2/2);
%     % caculat critical srain energy release rate for individual carbon
%     % fibres
% %     Wtotal_c(1:size(MminOverlap{i},3),i)=sum(sum(Wdeb_c+Wpo_c))/sum(sum(area_c));
%
%
%
% %     pause
%
%     %repeat for hybrid interfaces
%     diam_cg=bsxfun(@times,MminOverlapType{i}==1,(Dfg+Dfc));
%     area_cg=bsxfun(@times,MminOverlapType{i}==1,((pi*(Dfg/2)^2/2+pi*(Dfc/2)^2/2)/2/Vf));
%     Wdeb_cg=bsxfun(@times,pi.*diam_cg/4*Gm,MminOverlap{i});
%     Wpo_cg=bsxfun(@times,pi.*diam_cg/4*Tu,MminOverlap{i}.^2/2);
% %     Wtotal_cg(1:size(MminOverlap{i},3),i)=sum(sum(Wdeb_cg+Wpo_cg))/sum(sum(area_cg));
%
% %     pause
%
%     % repeat for glass interfaces
%     diam_g=bsxfun(@times,MminOverlapType{i}==0,Dfg);
%     area_g=bsxfun(@times,MminOverlapType{i}==0,(pi*(Dfg/2)^2/2/Vf));
%     Wdeb_g=bsxfun(@times,pi.*diam_g/4*Gm,MminOverlap{i});
%     Wpo_g=bsxfun(@times,pi.*diam_g/4*Tu,MminOverlap{i}.^2/2);
% %     Wtotal_g(1:size(MminOverlap{i},3),i)=sum(sum(Wdeb_g+Wpo_g))/sum(sum(area_g));
%
%
%
%     Wtotal_c(1:size(MminOverlap{i},3),i)=(sum(sum(Wdeb_c+Wpo_c))+(sum(sum(area_c))/(sum(sum(area_c))+sum(sum(area_g))))*sum(sum(Wdeb_cg+Wpo_cg)))/(sum(sum(area_c))+0.5*sum(sum(area_cg)));
%     Wtotal_g(1:size(MminOverlap{i},3),i)=(sum(sum(Wdeb_g+Wpo_g))+(sum(sum(area_g))/(sum(sum(area_c))+sum(sum(area_g))))*sum(sum(Wdeb_cg+Wpo_cg)))/(sum(sum(area_g))+0.5*sum(sum(area_cg)));
% % %  pause
% %
% %     Gcrit_c_cracks=(sum(sum(Wtotal_c))/sum(sum(area_c)))+0.5*(sum(sum(Wtotal_cg))/sum(sum(area_cg)));
% %     Gcrit_g(i)=(sum(sum(Wtotal_g))/sum(sum(area_g)))+0.5*(sum(sum(Wtotal_cg))/sum(sum(area_cg)));
%
% %     pause
% end
% Gcrit_homogenised=min(Wtotal);

%% method 3 for fracture toughness - joel's new analytical solution
% may need to change this if finding Gcrit from remaining overlap lengths
% (will be different for different RVEs)
Gcrit_c=findGcrit(lf,Vf,GiicmAvg,Tu0,dc,alphac,SmAvg,lcWei,XcWei,mcWei,0,0);
Gcrit_g=findGcrit(lf,Vf,GiicmAvg,Tu0,dg,alphag,SmAvg,lgWei,XgWei,mgWei,0,0);

%% record specimen properties
Specimen_dataset.Material_properties = Material_properties(Carbon_fibre,dc,ec,EcAvg,EcCoV,lcWei,mcWei,Sc,Glass_fibre,dg,eg,EgAvg,EgCoV,lgWei,mgWei,Sg,G,GiicmAvg,gmAvg,SmAvg,SmCoV,Tu0,Vc,Vf,lf,lspecimen,...
    pc,pcFailed,pg,pgFailed,pInterface,misalignment_switch,Gcrit_c,Gcrit_g);

% Numerical

if Fracture_method==1
    
    % Joel's version of fracture failure criterion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%FAILURE CRITERION %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dbstop if error
    % Numerical
    naGuess0=2; %the first value of na that will be tried
    
    for i=1:nRVE
        
        %EE{i},MSScurves{i},alphag(i),alphac(i)
        Xfibres=MSScurves{i}*Vf;
        Xmax1=max(SigSeries);
        Xmax2=mean(mean(max(Xfibres,[],3)));
        Xmax=mean(mean(max(Xfibres(2:end-1,2:end-1,:),[],3)));
        Xmax2/Xmax1
        Xmax/Xmax1
        
        iXmax=0;
        
        % Gcrit=Wtotal(end,i)
        
        Gcrit=Vc*Gcrit_c + (1-Vc)*Gcrit_g;
        
        % Calculate Jpp for the first guess
        %         [XcritGuess,JppGuess,imaxXRVEi,jmaxXRVEi]=SFC_fracture_Joel_hybrid_published(naGuess0,n_rows,n_columns,Xfibres(2:end-1,2:end-1,:),Df,Vf,Xmax,iXmax,Enom);
        [XcritGuess,JppGuess,imaxXRVEi,jmaxXRVEi]=SFC_fracture_Joel_hybrid_published(naGuess0,n_rows-2,n_columns-2,Xfibres(3:end-2,3:end-2,:),Df,Vf,Xmax,iXmax,Enom);
        isGuessCrit=JppGuess>Gcrit;
        
        % Calculate Jpp for the second guess
        naNow=naGuess0-isGuessCrit*2+1;
        [XcritNow,JppNow,imaxXRVEi,jmaxXRVEi]=SFC_fracture_Joel_hybrid_published(naNow,n_rows-2,n_columns-2,Xfibres(3:end-2,3:end-2,:),Df,Vf,Xmax,iXmax,Enom);
        isNowCrit=JppNow>Gcrit;
        
        naGuess=naGuess0;
        
        while isGuessCrit==isNowCrit
            % Update Guesses
            naGuess=naNow; XcritGuess=XcritNow;
            % Calculate new estimate
            naNow=naGuess-isGuessCrit*2+1;
            [XcritNow,JppNow,imaxXRVEi,jmaxXRVEi]=SFC_fracture_Joel_hybrid_published(naNow,n_rows-2,n_columns-2,Xfibres(3:end-2,3:end-2,:),Df,Vf,Xmax,iXmax,Enom);
            imaxXRVE{i}=imaxXRVEi;
            jmaxXRVE{i}=jmaxXRVEi;
            isNowCrit=JppNow>Gcrit;
        end
        imaxXRVE{i}=imaxXRVEi;
        jmaxXRVE{i}=jmaxXRVEi;
        
        % Assigns correct critical solution
        na(i)=isGuessCrit*naGuess+isNowCrit*naNow;
        Xcrit(i)=isGuessCrit*XcritGuess+isNowCrit*XcritNow;
        Specimen_dataset.RVE{i}.Xcrit=Xcrit(i);
        
        cluster_loc_counter=1;
        for i_cluster=1:na(i)
            for j_cluster=1:na(i)
                Specimen_dataset.RVE{i}.Nearly_critical_cluster_loc(cluster_loc_counter,1)=imaxXRVE{i}+i_cluster-1;
                Specimen_dataset.RVE{i}.Nearly_critical_cluster_loc(cluster_loc_counter,2)=jmaxXRVE{i}+j_cluster-1;
                cluster_loc_counter=cluster_loc_counter+1;
            end
        end
    end
    
    [Specimen_dataset.Specimen_fracture_strength, Specimen_dataset.Critical_RVE_index] = min(Xcrit);
    Specimen_dataset.Critical_cluster_loc=[];
    cluster_loc_counter=1;
    for i=1:na(Specimen_dataset.Critical_RVE_index)
        for j=1:na(Specimen_dataset.Critical_RVE_index)
            Specimen_dataset.Critical_cluster_loc(cluster_loc_counter,1)=imaxXRVE{Specimen_dataset.Critical_RVE_index}+i-1;
            Specimen_dataset.Critical_cluster_loc(cluster_loc_counter,2)=jmaxXRVE{Specimen_dataset.Critical_RVE_index}+j-1;
            cluster_loc_counter=cluster_loc_counter+1;
        end
    end
    
    % FJ - create duplicate stress-strain curve and apply toughness criteria
    SigSeries2=SigSeries;
    SigSeries2(find(SigSeries<min(Xcrit),1,'last'):end)=0; % Set end to 0
    % if length(SigSeries2)~=length(SigSeries)
    %    SigSeries2=SigSeries2(1:find(SigSeries2(2:end)==0)); % Remove end
    % end
    SigSeries2=SigSeries2(1:find(SigSeries2(2:end)==0)); % Remove end
    EpsSeries2=EpsSeries(1:length(SigSeries2));
    
    %% Fin's version
elseif Fracture_method == 2
    
    
    isGuessCrit=0;
    
    for i=1:nRVE
        
        %     EE{i},MSScurves{i},alphag(i),alphac(i)
        Xfibres=MSScurves{i}*Vf;
        Xmax1=max(SigSeries);
        Xmax2=mean(mean(max(Xfibres,[],3)));
        Xmax=mean(mean(max(Xfibres(2:end-1,2:end-1,:),[],3)));
        Xmax2/Xmax1;
        Xmax/Xmax1;
        
        iXmax=0;
        
        % FJ - first guess for Xcrit and J integral
        dbstop if error
        
        %     CrackSize_Rows=1; %the first value of crack size that will be tried
        %     CrackSize_Columns=1;
        %     Cumulative_cluster_data=[];
        
        
        %         [Xfracture(i),RVE_fracture_strain_increment(i),failed_cluster{i},Specimen_dataset.RVE{i}]=SFCfracture_EpsCrit_VariableShape_VariableGcrit_Clusters(Vc,nc,Gcrit_c,Gcrit_g,n_rows,n_columns,Xfibres(2:end-1,2:end-1,:),Df,dc,dg,Vf,Xmax,iXmax,Enom,Specimen_dataset.RVE{i});
        [Xfracture(i),RVE_fracture_strain_increment(i),failed_cluster{i},Specimen_dataset.RVE{i}]=SFCfracture_Clusters_EffectsOfRandomness_bwconncomp_othernewX(deltaEpsilon,EcAvg,EgAvg,Vc,nc,Gcrit_c,Gcrit_g,n_rows,n_columns,Xfibres(2:end-1,2:end-1,:),Df,dc,dg,Vf,Xmax,iXmax,Enom,Specimen_dataset.RVE{i},tm,lf,GiicmAvg,SmAvg,gmAvg,GiicmAvg,SmAvg,gmAvg,GiicmAvg,SmAvg,gmAvg);
        %     CutOffStress=50;
        %     [Xfracture(i),RVE_fracture_strain_increment(i),failed_cluster{i},Specimen_dataset.RVE{i}]=SFCfracture_Clusters_ReducedFibreArea_ResidualClusterStresses(Vc,nc,Gcrit_c,Gcrit_g,n_rows,n_columns,Xfibres(2:end-1,2:end-1,:),Df,dc,dg,Vf,Xmax,iXmax,Enom,Specimen_dataset.RVE{i},CutOffStress);
        
        
        Specimen_dataset.RVE{i}.Xcrit=Xfracture(i);
        Specimen_dataset.RVE{i}.Epscrit=RVE_fracture_strain_increment(i)*deltaEpsilon;
        Specimen_dataset.RVE{i}.Failed_cluster=failed_cluster{i};
        %         [Specimen_dataset.RVE{i}.Nearly_critical_cluster_loc(:,1),Specimen_dataset.RVE{i}.Nearly_critical_cluster_loc(:,2)]=ind2sub([size(Specimen_dataset.RVE{i}.RVE_fibre_type,1) size(Specimen_dataset.RVE{i}.RVE_fibre_type,2)],find(failed_cluster{i}==1));
        %     Specimen_dataset.RVE{i}.Reserve_factor_record=reserve_factor_record{i};
        %     Specimen_dataset.RVE{i}.Cluster_size_record=cluster_size_record{i};
        
        
        %     Eps_crit(i)=isGuessCrit*Eps_crit_guess+isNowCrit*Eps_crit_now;
        %     Specimen_dataset.RVE{i}.Xcrit=isGuessCrit*XcritGuess+isNowCrit*XcritNow;
        %     Specimen_dataset.RVE{i}.Gcrit=isGuessCrit*GcritGuess(i)+isNowCrit*GcritNow(i);
        %     Specimen_dataset.RVE{i}.Cluster_data=isGuessCrit*Cluster_data_guess+isNowCrit*Cluster_data_now;
        %     Specimen_dataset.RVE{i}.Cumulative_cluster_data=Cumulative_cluster_data;
        %     [Fibre_failure_increments Fibre_failure_increments]=max(Xfibres(2:end-1,2:end-1,:),[],3);
        %     Specimen_dataset.RVE{i}.Fibre_failure_increments=Fibre_failure_increments;
        
        
        % Don't save fibre stresses for each RVE unless debugging - it takes up
        % soooo much RAM!
        %         Specimen_dataset.RVE{i}.Fibre_stresses=Xfibres(2:end-1,2:end-1,:);
        
        % for interpolation 5 method - if load drop and interpolating on
        % strain, then fracture will occur at higher stress
        
        if interpolation_method ==5
            [RVE_max_stress(i) RVE_max_stress_increment(i)] = max(Specimen_dataset.RVE{i}.Stress_strain_unprocessed(1:RVE_fracture_strain_increment(i),2));
            % FJ - change the fracture condition flag and set the critical
            % cluster to the full cross-section & also  cluster size?
            if  RVE_max_stress(i) > Xfracture(i)
                Xfracture(i)=RVE_max_stress(i);
                RVE_fracture_strain_increment(i)=RVE_max_stress_increment(i);
                
            elseif RVE_max_stress(i) == Xfracture(i) && RVE_max_stress_increment(i) < RVE_fracture_strain_increment(i)
                Xfracture(i)=RVE_max_stress(i);
                RVE_fracture_strain_increment(i)=RVE_max_stress_increment(i);
                
            end
        end
    end
    
    if interpolation_method ==5
        [Specimen_dataset.Specimen_fracture_strength, Specimen_dataset.Critical_RVE_index] = min(Xfracture);
        % FJ - this method doesn't work if large level of misalignment
        Specimen_fracture_strain_increment = find(SigSeries>=Specimen_dataset.Specimen_fracture_strength,1,'first');
        % sometimes it falls over due to the slight inaccuracies from the
        % interpolation
        if isempty(Specimen_fracture_strain_increment)==1
            Specimen_fracture_strain_increment=size(SigSeries,1);
        end
    else
        [Specimen_fracture_strain_increment critical_RVE_index]=min(RVE_fracture_strain_increment);
        Specimen_dataset.Specimen_fracture_strength=Xfracture(critical_RVE_index);
        Specimen_dataset.Critical_RVE_index=critical_RVE_index;
    end
    Specimen_dataset.Specimen_fracture_strain=(Specimen_fracture_strain_increment-1)*deltaEpsilon;
    Specimen_dataset.Critical_cluster=failed_cluster{Specimen_dataset.Critical_RVE_index};
    if isempty(Specimen_dataset.Critical_cluster)~=1
        [Specimen_dataset.Critical_cluster_loc(:,1),Specimen_dataset.Critical_cluster_loc(:,2)]=ind2sub([size(Specimen_dataset.RVE{Specimen_dataset.Critical_RVE_index}.RVE_fibre_type,1) size(Specimen_dataset.RVE{Specimen_dataset.Critical_RVE_index}.RVE_fibre_type,2)],find(failed_cluster{Specimen_dataset.Critical_RVE_index}==1));
        Specimen_dataset.Critical_cluster_centroid(1,1)=round(sum(Specimen_dataset.Critical_cluster_loc(:,1),1)/size(Specimen_dataset.Critical_cluster_loc(:,1),1)+rand()-0.5);
        Specimen_dataset.Critical_cluster_centroid(1,2)=round(sum(Specimen_dataset.Critical_cluster_loc(:,2),1)/size(Specimen_dataset.Critical_cluster_loc(:,2),1)+rand()-0.5);
        Specimen_dataset.Crit_damage_threshold=Specimen_dataset.RVE{Specimen_dataset.Critical_RVE_index}.Crit_damage_threshold;
        Specimen_dataset.Crit_cluster_size=Specimen_dataset.RVE{Specimen_dataset.Critical_RVE_index}.Cluster_size_record;
    end
    
    
    % FJ - plot results
    % figure()
    % hold on
    % plot(EpsSeries,SigSeries)
    % xlabel('Strain (\%)','interpreter','latex')
    % ylabel('Stress (MPa)','interpreter','latex')
    % set(gca,'TickLabelInterpreter', 'latex');
    % legend boxoff
    % xlim([0 1])
    % ylim([0 750])
    
    %% althernativemethods to combine RVEs
    % FJ - find the max strength of the strongest RVE
    deltaSigma_NewFracture=3;
    
    maxi=0;
    for i=1:nRVE
        maxi=max(maxi,max(EE{i}(:,2)));
    end
    maxi;
    
    % % FJ - find the strength of the weakest RVE
    % mini=maxi;
    % for i=1:nRVE
    %     mini=min(mini,max(EE{i}(:,2)));
    % end
    % mini;
    
    % FJ - set strength scale using highest strength value from all RVEs
    % FJ - check whether this should be based on max RVE strength or min RVE
    % strength. Why would it be different between stress and strain?
    NewSigmaGrid_NewFracture=[0:deltaSigma_NewFracture:maxi+deltaSigma_NewFracture]';
    
    
    % FJ - calculate number of strain points by dividing the weakest RVE by the
    % strength increment and then round down to closest integer
    % FJ - EpsSeries is quite sensitive to strain, as gstiffness is low in this region, the +1 at the end makes a large difference! Also this is a global stress-strin relation, try find on x and y
    % neps=find(EpsSeries>(Specimen_fracture_strain_increment*deltaEpsilon),1,'first')+1;
    
    %     neps=floor(mini/deltaSigma_NewFracture);
    % FJ - pre-allocate Eps outside of loop
    Eps=[];
    
    % FJ - interpolate strains to give RVE stress-strain response up to
    % lowest RE strength with common stress-strain scale across RVEs
    for i=1:nRVE
        % FJ - define x and y as a vector of RVE strains and RVE strengths respectively - this is overwritten with
        % each loop
        x=EE{i}(:,1);
        y=EE{i}(:,2);
        
        
        % FJ - trim the contents of the matrices x and y up until the strain at
        % which the weakst RVE is at its max
        x=x(1:neps);
        y=y(1:neps);
        
        
        
        x2=x;
        y2=y;
        
        % Make x and y monotonically increasing after fracture strain, in order
        % to ensure a stress plateau at the correct height - doesn't work with
        % cases where load drop is recovered!
        %     x(Specimen_fracture_strain_increment+1:end);
        %     y(Specimen_fracture_strain_increment+1:end)=[(y(Specimen_fracture_strain_increment)+0.01):0.01:((length(y)-(Specimen_fracture_strain_increment))*0.01+y(Specimen_fracture_strain_increment))];
        
        %Smaller epsilon
        n=0;
        try
            while n<length(y)-1
                % TESTING
                %             if y(end-(n+1))>y(end-n) && x(end-n)<=((Specimen_fracture_strain_increment-10)*deltaEpsilon)
                if y(end-(n+1))>y(end-n) && x(end-n)<=((Specimen_fracture_strain_increment)*deltaEpsilon)
                    %             if y(end-(n+1))>y(end-n)
                    y=[y(1:end-(n+1)-1); y(end-(n+1)+1:end)];
                    x=[x(1:end-(n+1)-1); x(end-(n+1)+1:end)];
                else
                    n=n+1;
                end
            end
        catch e
            %             save(['VTF_error_Vc_' num2str(Vc) '_' filename],'Specimen_dataset','e','y2','y','x2','x','n','Specimen_fracture_strain_increment','deltaEpsilon','-v7.3')
            save(['VTF_' defect_case '_error_Vc_' num2str(Vc) '_' filename],'-v7.3')
            error(['error in file: ' e.stack(1).name ', line: ' num2str(e.stack(1).line) ': ' e.message]);
        end
        
        x3=x;
        y3=y;
        
        try
            % now make monotonically increasing for rest of curve (bottom to
            % top)
            % TESTING
            %         CutOffIncrement=find(x>=(Specimen_fracture_strain_increment-10+1)*deltaEpsilon,1,'first')
            CutOffIncrement=find(x>=(Specimen_fracture_strain_increment+1)*deltaEpsilon,1,'first');
            y(CutOffIncrement:end)=y(CutOffIncrement-1)+[deltaSigma_NewFracture:deltaSigma_NewFracture:(length(y)-(CutOffIncrement-1))*deltaSigma_NewFracture];
        catch e
            %             save(['VTF_error_Vc_' num2str(Vc) '_' filename],'Specimen_dataset','e','y3','y','x3','x','CutOffIncrement','Specimen_fracture_strain_increment','deltaEpsilon','deltaSigma_NewFracture','-v7.3')
            save(['VTF_' defect_case '_error_Vc_' num2str(Vc) '_' filename],'-v7.3')
            error(['error in file: ' e.stack(1).name ', line: ' num2str(e.stack(1).line) ': ' e.message]);
        end
        
        %         n=1;
        %         while n<length(x)-1
        %             if y(n+1)<y(n)
        %                 y=[y(1:n); y(n+2:end)];
        %                 x=[x(1:n); x(n+2:end)];
        %             else
        %                 n=n+1;
        %             end
        %         end
        
        
        %     figure()
        %     hold on
        %     plot(x,y)
        %     xlabel('Strain (\%)','interpreter','latex')
        %     ylabel('Stress (MPa)','interpreter','latex')
        %     set(gca,'TickLabelInterpreter', 'latex');
        %     legend boxoff
        %     xlim([0 1])
        %     ylim([0 1000])
        
        % FJ - interpolate data using common scales for all RVEs
        %     NewEpsilonGrid = x;
        if max(isnan(x))==1 || max(isnan(y))==1 || max(isnan(NewSigmaGrid_NewFracture))==1
            1;
        end
        NewEpsilonGrid_NewFracture = interp1(y,x,NewSigmaGrid_NewFracture);
        %         NewEpsilonGrid = interpn(y,x,NewSigmaGrid,'cubic');
        % FJ - trim the strain vector to the number of strain points
        %     NewEpsilonGrid_NewFracture=NewEpsilonGrid_NewFracture(1:find(EpsSeries>(Specimen_fracture_strain_increment*deltaEpsilon),1,'first')+1);
        
        %     Strengths(:,i)=y(1:neps);
        
        % FJ - populate matrix of strains
        Eps_NewFracture(:,i)=NewEpsilonGrid_NewFracture;
    end
    
    
    % FJ - calculate the average strain vector across all RVEs
    EpsSeries_NewFracture=sum(Eps_NewFracture,2)/nRVE;
    % EpsSeries_NewFracture=EpsSeries_NewFracture(1:(find(EpsSeries>(Specimen_fracture_strain_increment*deltaEpsilon),1,'first')+1));
    EpsSeries_NewFracture=EpsSeries_NewFracture(1:(find(EpsSeries_NewFracture>((Specimen_fracture_strain_increment-1)*deltaEpsilon),1,'first')));
    % EpsSeries=x(1:neps);
    % FJ - trim the strength vector so that the max strength is that of the
    % weakest RVE (chain theory)
    SigSeries_NewFracture=NewSigmaGrid_NewFracture(1:length(EpsSeries_NewFracture));
    % SigSeries=sum(Strengths,2)/nRVE;
    
    
    
    if interpolation_method == 5
        % FJ - create duplicate stress-strain curve and apply toughness criteria
        SigSeries2=SigSeries;
        SigSeries2=SigSeries2(1:(find(SigSeries2>=Specimen_dataset.Specimen_fracture_strength-1,1,'first')));
        % FJ - trim the strength vector so that the max strength is that of the
        % weakest RVE (chain theory)
        EpsSeries2=EpsSeries(1:length(SigSeries2));
        
    else
        % FJ - create duplicate stress-strain curve and apply toughness criteria
        EpsSeries2=EpsSeries;
        EpsSeries2=EpsSeries2(1:(find(EpsSeries2>=((Specimen_fracture_strain_increment-1)*deltaEpsilon),1,'first')));
        % FJ - trim the strength vector so that the max strength is that of the
        % weakest RVE (chain theory)
        SigSeries2=SigSeries(1:length(EpsSeries2));
        
        % EpsSeries2(find(EpsSeries>((Specimen_fracture_strain_increment-1)*deltaEpsilon),1,'first'):end)=0; % Set end to 0
        % if length(EpsSeries2)~=length(EpsSeries)
        %     EpsSeries2=EpsSeries2(1:find(EpsSeries2(2:end)==0)); % Remove end
        % end
        % % EpsSeries2=EpsSeries2(1:find(EpsSeries2(2:end)==0)); % Remove end
        %
        % EpsSeries2=EpsSeries(1:(find(EpsSeries>(Specimen_fracture_strain_increment*deltaEpsilon),1,'first')+1));
        %
        % SigSeries2=SigSeries(1:length(EpsSeries2));
        
        
        % EpsSeries2=EpsSeries;
        % EpsSeries2=EpsSeries2(1:Specimen_fracture_strain_increment);
        % SigSeries2=SigSeries(1:Specimen_fracture_strain_increment);
        
    end
    
    %% comment this out if using joels fracture criterion
    for i=1:nRVE
        Specimen_dataset.RVE{i}.Stress_strain_fracture_unprocessed=Specimen_dataset.RVE{i}.Stress_strain_unprocessed(1:RVE_fracture_strain_increment(i),:);
        
    end
end



% Store stress-strain response for plotting
Specimen_SS_Response=[EpsSeries,SigSeries];
Specimen_dataset.Stress_strain=Specimen_SS_Response;

if interpolation_method == 1
    % Specimen_SS_Response2=[EpsSeries_NewFracture,SigSeries_NewFracture];
    Specimen_SS_Response2=[EpsSeries_NewFracture,SigSeries_NewFracture];
    Specimen_dataset.Stress_strain_fracture=Specimen_SS_Response2;
else
    Specimen_SS_Response2=[EpsSeries2,SigSeries2];
    Specimen_dataset.Stress_strain_fracture=Specimen_SS_Response2;
end

if strcmp(Alignment_method,'None')~=1
    Constituent_props=Constituent_properties(6.0,6.0,G/1000,0.45,EgAvg/1000,6.0,5.1,0.28,EcAvg/1000,6,5.1,0.28,Vf,Vc);
    
%     % laimina properties taken from Joel's hybrid paper
%     switch Composite_type
%         case 'HSC'
%             Constituent_props=Constituent_properties(4.0,4.0,G/1000,0.35,EgAvg/1000,EgAvg/1000,30,0.2,EcAvg/1000,20,27,0.2,Vf,Vc);
%         case 'HMC_EG'
%             Constituent_props=Constituent_properties(4.0,4.0,G/1000,0.35,EgAvg/1000,EgAvg/1000,30,0.2,EcAvg/1000,7,32,0.2,Vf,Vc);
%         otherwise
%             error('sorry, lamina properties are not defined for other material types');
%     end
    
    [Specimen_dataset]= EffectsOfRandomness_MisAlignment(Specimen_dataset,Constituent_props,Alignment_distribution,Sigma_alignment,Alignment_method);
    Specimen_SS_Response2=Specimen_dataset.Misaligned_stress_strain;
    
end

% store specimen preformance chics
Specimen_dataset.Gcrit=Specimen_dataset.RVE{Specimen_dataset.Critical_RVE_index}.Gcrit;
Specimen_dataset.Initial_stiffness=(Specimen_SS_Response2(2,2)/Specimen_SS_Response2(2,1))*100;
Specimen_dataset.Ultimate_strain=Specimen_SS_Response2(end,1);
Specimen_dataset.Ultimate_strength=Specimen_SS_Response2(end,2);

switch failure_modes_switch
    case 1
        Specimen_dataset=determine_fibre_failure_modes(Specimen_dataset,deltaEpsilon);
    otherwise
end

Strain_offset_yield_projection=((Specimen_SS_Response2(:,1)-0.1)./100).*Specimen_dataset.Initial_stiffness;

Specimen_dataset.Yield_strength=Specimen_SS_Response2(find(Strain_offset_yield_projection>=Specimen_SS_Response2(:,2),1,'first'),2);

Specimen_dataset.Pseudo_ductile_strain=Specimen_dataset.Ultimate_strain-(Specimen_dataset.Ultimate_strength/Specimen_dataset.Initial_stiffness*100);

Specimen_dataset.Max_load_drop=determine_max_load_drop(Specimen_dataset.Stress_strain_fracture);

save(['VTF_' single_multi_flag '_' objective_fun_flag '_' defect_case '_run_' filename],'Specimen_dataset','-v7.3')

disp('specimen complete')

% write dummy file to workspace to tell HPC when it's successfully finished
fileID = fopen('sentinel.txt','w');
nbytes = fprintf(fileID,'specimen_complete');
end

function Specimen_dataset=determine_fibre_failure_modes(Specimen_dataset,deltaEpsilon)
Critical_RVE_index=Specimen_dataset.Critical_RVE_index;
Fragmentation_failures=zeros(size(Specimen_dataset.RVE{Critical_RVE_index}.Failure_mode_record));
Debonding_failures=Fragmentation_failures;
Softening_failures=Fragmentation_failures;
Fibre_failures=Fragmentation_failures;
for ii=1:size(Specimen_dataset.RVE{Critical_RVE_index}.Failure_mode_record,1)
    for jj=1:size(Specimen_dataset.RVE{Critical_RVE_index}.Failure_mode_record,2)
        for kk=1:length(find(Specimen_dataset.RVE{Critical_RVE_index}.Failure_mode_record{ii,jj}{1}(:,1)<=(Specimen_dataset.Specimen_fracture_strain/deltaEpsilon)))
            Fragmentation_failures(ii,jj)=Fragmentation_failures(ii,jj)+(Specimen_dataset.RVE{Critical_RVE_index}.Failure_mode_record{ii,jj}{2}{kk}==2);
            Debonding_failures(ii,jj)=Debonding_failures(ii,jj)+(Specimen_dataset.RVE{Critical_RVE_index}.Failure_mode_record{ii,jj}{2}{kk}==3);
            Softening_failures(ii,jj)=Softening_failures(ii,jj)+(Specimen_dataset.RVE{Critical_RVE_index}.Failure_mode_record{ii,jj}{2}{kk}==4);
            Fibre_failures(ii,jj)=Fibre_failures(ii,jj)+(Specimen_dataset.RVE{Critical_RVE_index}.Failure_mode_record{ii,jj}{2}{kk}==1);
        end
    end
end

Specimen_dataset.Fragmentation_failures=Fragmentation_failures;
Specimen_dataset.Debonding_failures=Debonding_failures;
Specimen_dataset.Softening_failures=Softening_failures;
Specimen_dataset.Fibre_failures=Fibre_failures;

end

function Max_load_drop = determine_max_load_drop(Stress_strain)
% finds the maximum load drop in the final stress-strain curve (part of
% optimisation objective function)
[~,~,~,load_drops]=findpeaks(Stress_strain(:,2));

if isempty(load_drops)==1 || max(load_drops)<=0
    Max_load_drop=0;
else
    Max_load_drop=max(load_drops);
end

end