function results = Initialise_VTF_BayesianOptimizationObject(single_multi_flag,objective_fun_flag,batch_index,num_runs,exploration_ratio,pc,pInterface,pcFailed,pg,pgFailed,misalignment_switch,rng_switch)

% use batch index to describe which stress-strain curves features to
% optimise
exploration_ratio=exploration_ratio/100;

if max(abs([pc,pInterface,pcFailed,pg,pgFailed,misalignment_switch]))>0
    defect_case='defective';
else
    defect_case='pristine';
end

switch objective_fun_flag
    case 'single'
        wInitS=batch_index==1;
        wUltStrain=batch_index==2;
        wPdStrain=batch_index==3;
        wUltStrength=batch_index==4;
        wYStrength=batch_index==5;
    case 'multiplicative'
        wInitS=batch_index==1;
        wUltStrain=batch_index==3;
        wPdStrain=batch_index==2;
        wUltStrength=batch_index==1 || batch_index==2;
        wYStrength=batch_index==3;
    case 'additive'     
        wInitS=(batch_index>=1 && batch_index<=101).*(1-(batch_index-1)./100);
        wUltStrength=(batch_index>=1 && batch_index<=101).*((batch_index-1)./100)+(batch_index>=102 && batch_index<=201).*(2-(batch_index-1)./100);
        wPdStrain=(batch_index>101 && batch_index<=201).*((batch_index-1)./100-1);
        wYStrength=(batch_index>=202 && batch_index<=302).*(3-(batch_index-2)./100);
        wUltStrain=(batch_index>=202 && batch_index<=302).*((batch_index-2)./100-2);
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
[initial_data,size_of_data]=determine_initial_values();

% decalre handle to objective function
fun = @(x)VTF_objective_function(objective_fun_flag,x.Carbon_fibre,x.Glass_fibre,x.Vc,x.lf,x.SmAvg,x.G,x.GiicmAvg,wInitS,wUltStrain,wPdStrain,wUltStrength,wYStrength,pc,pInterface,pcFailed,pg,pgFailed,misalignment_switch,rng_switch,[],[],[],[],[]);


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
    'PlotFcn',{@plotMinObjective,@plotConstraintModels,@plotObjectiveModel});

% save results
save([single_multi_flag '_' objective_fun_flag '_' defect_case '_results_' num2str(wInitS*100) '_' num2str(wUltStrain*100) '_' num2str(wPdStrain*100) '_' num2str(wUltStrength*100) '_' num2str(wYStrength*100) '.mat'],'results','-v7.3');

end

% determine initial values of inputs and objective function outputs from
% previous datasets
function [initial_data,size_of_data]=determine_initial_values()

Carbon_fibre=categorical({'T1000GB','K13D','XN-05','T300','M40B','M60JB','XN-90','P120J','P75S','C320','T800H','HTA5131','C124','C100','FliteStrand_S_ZT','GF'});
Glass_fibre=categorical({'T1000GB','K13D','XN-05','T300','M40B','M60JB','XN-90','P120J','P75S','C320','T800H','HTA5131','C124','C100','FliteStrand_S_ZT','GF'});
size_of_data=size(Carbon_fibre,2); % record size of dataset
Vc=rand(size_of_data,1);
lf=(rand(size_of_data,1)*11.5+0.5)*1000; % convert lf to um - might change this in the future
SmAvg=rand(size_of_data,1)*60+40;
G=rand(size_of_data,1)*800+1000;
GiicmAvg=rand(size_of_data,1)*0.4+0.6;

initial_data=table(Carbon_fibre',flipud(Glass_fibre'),Vc,lf,SmAvg,G,GiicmAvg); % put inputs into a table

size_of_data=size(Vc,1); % record size of dataset

end