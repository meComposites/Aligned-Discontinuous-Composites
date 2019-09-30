function results = Create_VTF_BayesianOptimizationObject(single_multi_flag,objective_fun_flag,objective_function_1,num_runs,exploration_ratio,pc,pInterface,pcFailed,pg,pgFailed,misalignment_switch,rng_switch)

% use batch index to describe which stress-strain curves features to
% optimise
exploration_ratio=exploration_ratio/100;

if max(abs([pc,pInterface,pcFailed,pg,pgFailed,misalignment_switch]))>0
    defect_case='defective';
else
    defect_case='pristine';
end

[wInitS,wUltStrain,wPdStrain,wUltStrength,wYStrength]=deal(0);
switch objective_fun_flag
    case 'single'
        switch objective_function_1
            case{'Initial_stiffness'}
                wInitS=1;
            case{'Ultimate_strain'}
                wUltStrain=1;
            case{'Pseudo_ductile_strain'}
                wPdStrain=1;
            case{'Ultimate_strength'}
                wUltStrength=1;
            case{'Yield_strength'}
                wYStrength=1;
            otherwise
                error("incorrect objective function 1 name. Please choose from 'Initial_stiffness', 'Ultimate_strain', 'Pseudo_ductile_strain', 'Ultimate_strength', or 'Yield_strength'.");
        end
           
%         wInitS=batch_index==1;
%         wUltStrain=batch_index==2;
%         wPdStrain=batch_index==3;
%         wUltStrength=batch_index==4;
%         wYStrength=batch_index==5;
%     case 'multiplicative'
%         wInitS=batch_index==1;
%         wUltStrain=batch_index==2;
%         wPdStrain=batch_index==3;
%         wUltStrength=batch_index==1 || batch_index==2;
%         wYStrength=batch_index==3;
%     case 'additive'     
%         wInitS=(batch_index>=1 && batch_index<=101).*(1-(batch_index-1)./100);
%         wUltStrength=(batch_index>=1 && batch_index<=101).*((batch_index-1)./100)+(batch_index>=102 && batch_index<=201).*(2-(batch_index-1)./100);
%         wUltStrain=(batch_index>101 && batch_index<=201).*((batch_index-1)./100-1);
%         wYStrength=(batch_index>=202 && batch_index<=302).*(3-(batch_index-2)./100);
%         wPdStrain=(batch_index>=202 && batch_index<=302).*((batch_index-2)./100-2);
%     case 'additive_detailed'
%         wInitS=0;
%         wUltStrength=(batch_index>=1 && batch_index<=201).*(0.75-(batch_index-1)./400);
%         wUltStrain=(batch_index>=1 && batch_index<=201).*(0.25+(batch_index-1)./400);
%         wYStrength=(batch_index>=202 && batch_index<=402).*(1.25-(batch_index-2)./400);
%         wPdStrain=(batch_index>=202 && batch_index<=402).*(-0.25+(batch_index-2)./400);
%         if batch_index >=1 && batch_index <= 201
%             objective_fun_flag='additive_detailed_ultStrengthUltStrain';
%         else
%             objective_fun_flag='additive_detailed_yStrengthPdStrain';
%         end
    otherwise
        error('incorrect objective function selection')
end

% set rng generator to random SIMD twister
switch rng_switch
    case 1
        rng('shuffle','simdTwister')
    otherwise
        rng(rng_switch)
end

% set optimizable variables
% With Carbon_fibre and Glass Fibre defined as such, you will produce "resultsNoSort.mat"
% Carbon_fibre = optimizableVariable('Carbon_fibre',{'T1000GB','K13D','XN-05','T300','M40B','M60JB','XN-90','P120J','P75S','C320','T800H','HTA5131','C124','C100','FliteStrand_S_ZT','GF'},'Type','categorical');
% Glass_fibre = optimizableVariable('Glass_fibre',{'T1000GB','K13D','XN-05','T300','M40B','M60JB','XN-90','P120J','P75S','C320','T800H','HTA5131','C124','C100','FliteStrand_S_ZT','GF'},'Type','categorical');
% Use this alternative definition of Carbon_fibre and Glass_fibre to create
% "resultsWithSort.mat"
Carbon_fibre = optimizableVariable('Carbon_fibre',sort({'T1000GB','K13D','XN-05','T300','M40B','M60JB','XN-90','P120J','P75S','C320','T800H','HTA5131','C124','C100','FliteStrand_S_ZT','GF'}),'Type','categorical');
Glass_fibre = optimizableVariable('Glass_fibre',sort({'T1000GB','K13D','XN-05','T300','M40B','M60JB','XN-90','P120J','P75S','C320','T800H','HTA5131','C124','C100','FliteStrand_S_ZT','GF'}),'Type','categorical');

% definition of remaining optimizable variables
Vc = optimizableVariable('Vc',[0,1],'Type','real');
lf = optimizableVariable('lf',[500,12000],'Type','real');
% Matrix = optimizableVariable('Matrix',{'Epoxy','Polypropylene','PEEK'},'Type','categorical');
SmAvg = optimizableVariable('SmAvg',[40,100],'Type','real');
G = optimizableVariable('G',[1000,1800],'Type','real');
GiicmAvg = optimizableVariable('GiicmAvg',[0.6,1.0],'Type','real');

% definie initial data for optimisation according to what we want to
% optimise - varargout from function may leave values undefined - make sure
% objectve function eneds with min values (for additive MOBO)
[initial_data,initial_objective,initial_constraints,initial_errors,size_of_data,minInitS,minPdStrain,minUltStrain,minUltStrength,minYStrength]=determine_initial_values(defect_case,objective_fun_flag,wInitS,wUltStrain,wPdStrain,wUltStrength,wYStrength);

% decalre handle to objective function
fun = @(x)VTF_objective_function(single_multi_flag,objective_fun_flag,x.Carbon_fibre,x.Glass_fibre,x.Vc,x.lf,x.SmAvg,x.G,x.GiicmAvg,wInitS,wUltStrain,wPdStrain,wUltStrength,wYStrength,pc,pInterface,pcFailed,pg,pgFailed,misalignment_switch,rng_switch,minInitS,minPdStrain,minUltStrain,minUltStrength,minYStrength);


% run bayesian optimisation with the following settings
results = bayesopt(fun,[Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg],...
    'Verbose',2,...
    'AcquisitionFunctionName','expected-improvement-plus',...
    'IsObjectiveDeterministic',0,...
    'ExplorationRatio',exploration_ratio,...
    'MaxObjectiveEvaluations',size_of_data+num_runs,... % run the optimisation for num runs = size of input data + additional required runs (set num_runs to zero to bypass future evaluations of objective function)
    'NumCoupledConstraints',1,...
    'AreCoupledConstraintsDeterministic',false,...% constraint itself is deterministic, but output from model is stochastic? - non-deterministic coupled constraint seems to work well
    'InitialX',initial_data,...% initial inputs from previous data
    'InitialObjective',initial_objective,...
    'InitialConstraintViolations',initial_constraints,...
    'InitialErrorValues',initial_errors,...
    'NumSeedPoints',1,...
    'PlotFcn',{@plotMinObjective,@plotConstraintModels,@plotObjectiveModel});

% save results
save([single_multi_flag '_' objective_fun_flag '_' defect_case '_results_' objective_function_1 '.mat'],'results','-v7.3');

end

% determine initial values of inputs and objective function outputs from
% previous datasets
function [initial_data,initial_objective,initial_constraints,initial_errors,size_of_data,minInitS,minPdStrain,minUltStrain,minUltStrength,minYStrength]=determine_initial_values(defect_case,obective_fun_flag,wInitS,wUltStrain,wPdStrain,wUltStrength,wYStrength)

try
%     load([defect_case '_VTF_data.mat'])
    load([defect_case '_VTF_data_corrected_test.mat'])
catch
    initial_data=[];initial_objective=[];initial_constraints=[];initial_errors=[];size_of_data=0;minInitS=[];minPdStrain=[];minUltStrain=[];minUltStrength=[];minYStrength=[];
    disp(['no observed VTF data - time to start from scratch!'])
    return
end

lf=lf*1000; % convert lf to um - might change this in the future
initial_data=table(Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg); % put inputs into a table
% determine value of yield strength
B=zeros(length(Vc),1);
B(isnan(Yield_strength)==1)=Ultimate_strength(isnan(Yield_strength)==1);
B(isnan(Yield_strength)==0)=Yield_strength(isnan(Yield_strength)==0);

% evaluate initial objective for each of the inputs, according to weights
% and type of objective function that's being analysed
switch obective_fun_flag
    case 'single'
        initial_objective=-1*(wInitS*(Initial_stiffness)+wUltStrain*(Ultimate_strain)+wPdStrain*(Pseudo_ductile_strain)+wUltStrength*(Ultimate_strength)+wYStrength*(B));
        minInitS=[];minPdStrain=[];minUltStrain=[];minUltStrength=[];minYStrength=[];
    case 'multiplicative'
        initial_objective=-1*(wInitS.*(Initial_stiffness).*wUltStrength.*(Ultimate_strength)+wUltStrength.*(Ultimate_strength).*wUltStrain.*(Ultimate_strain)+wYStrength.*(B).*wPdStrain.*(Pseudo_ductile_strain));
        minInitS=[];minPdStrain=[];minUltStrain=[];minUltStrength=[];minYStrength=[];
    case {'additive','additive_detailed_ultStrengthUltStrain','additive_detailed_yStrengthPdStrain'}
        switch defect_case
            case 'pristine'
                load('SOBO_FibreType_MinObesrved.mat')
            case 'defective'
                load('SOBO_defective_FibreType_MinObesrved.mat')
        end
        initial_objective=-1.*(wInitS.*(Initial_stiffness./minInitS)+wUltStrain.*(Ultimate_strain./minUltStrain)+wPdStrain.*(Pseudo_ductile_strain./minPdStrain)+wUltStrength.*(Ultimate_strength./minUltStrength)+wYStrength.*(B./minYStrength));
end

% set constraints and error traces
initial_constraints=-50+Max_load_drop;
initial_errors=2*isnan(initial_objective)-1;

size_of_data=size(Vc,1); % record size of dataset

end