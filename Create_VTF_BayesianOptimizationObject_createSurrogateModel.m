function results = Create_VTF_BayesianOptimizationObject_createSurrogateModel(e_constraint_name,e_constraint_min,n_e_constraint_slices,e_constraint_max,objective_name,single_multi_flag,objective_fun_flag,batch_index,num_runs,exploration_ratio,pc,pInterface,pcFailed,pg,pgFailed,misalignment_switch,rng_switch)

if batch_index>n_e_constraint_slices
    error('too many batch indices')
end

% use batch index to describe which stress-strain curves features to
% optimise
exploration_ratio=exploration_ratio/100;

if max(abs([pc,pInterface,pcFailed,pg,pgFailed,misalignment_switch]))>0
    defect_case='defective';
else
    defect_case='pristine';
end

switch objective_fun_flag
    case {'eConstraint','eConstraintHybrid'}
        wInitS=0;
        wUltStrain=0;
        wPdStrain=0;
        wUltStrength=0;
        wYStrength=0;
        eval([objective_name '=1;']);
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
[initial_data,initial_objective,initial_constraints,initial_errors,size_of_data,minInitS,minPdStrain,minUltStrain,minUltStrength,minYStrength]=determine_initial_values(e_constraint_name,e_constraint_min,n_e_constraint_slices,e_constraint_max,batch_index,defect_case,objective_fun_flag,wInitS,wUltStrain,wPdStrain,wUltStrength,wYStrength);

% decalre handle to objective function
fun = @(x)VTF_objective_function_eConstraint_hybrid(single_multi_flag,objective_fun_flag,objective_name,e_constraint_name,e_constraint_min,n_e_constraint_slices,e_constraint_max,batch_index,x.Carbon_fibre,x.Glass_fibre,x.Vc,x.lf,x.SmAvg,x.G,x.GiicmAvg,wInitS,wUltStrain,wPdStrain,wUltStrength,wYStrength,pc,pInterface,pcFailed,pg,pgFailed,misalignment_switch,rng_switch,minInitS,minPdStrain,minUltStrain,minUltStrength,minYStrength);
 
 % relax lower constraint if on first slice
if batch_index == 1
    numCoupledConstraints=1;
    areCoupledConstraintsDeterministic=[false];
else
    numCoupledConstraints=2;
    areCoupledConstraintsDeterministic=[false,false];
end

results = bayesopt(fun,[Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg],...
    'Verbose',2,...
    'AcquisitionFunctionName','expected-improvement-plus',...
    'IsObjectiveDeterministic',0,...
    'ExplorationRatio',exploration_ratio,...
    'MaxObjectiveEvaluations',size_of_data+num_runs,... % run the optimisation for num runs = size of input data + additional required runs (set num_runs to zero to bypass future evaluations of objective function)
    'NumCoupledConstraints',numCoupledConstraints,...
    'AreCoupledConstraintsDeterministic',areCoupledConstraintsDeterministic,...% constraint itself is deterministic, but output from model is stochastic? - non-deterministic coupled constraint seems to work well
    'InitialX',initial_data,...% initial inputs from previous data
    'InitialObjective',initial_objective,...
    'InitialConstraintViolations',initial_constraints,...
    'InitialErrorValues',initial_errors,...
    'NumSeedPoints',1,... % can this be set to zero?
    'PlotFcn',{@plotMinObjective,@plotConstraintModels,@plotObjectiveModel});    

% save results
save([single_multi_flag '_' objective_fun_flag '_' defect_case '_results_max_' objective_name '_contrained_' e_constraint_name '_segment_' num2str(batch_index) '.mat'],'results','-v7.3');

% write dummy file to workspace to tell HPC when it's successfully finished
fileID = fopen('sentinel.txt','w');
nbytes = fprintf(fileID,'optimisation complete');
end

% determine initial values of inputs and objective function outputs from
% previous datasets
function [initial_data,initial_objective,initial_constraints,initial_errors,size_of_data,minInitS,minPdStrain,minUltStrain,minUltStrength,minYStrength]=determine_initial_values(e_constraint_name,e_constraint_min,n_e_constraint_slices,e_constraint_max,batch_index,defect_case,obective_fun_flag,wInitS,wUltStrain,wPdStrain,wUltStrength,wYStrength)

% VTF_data_name=['complete_' defect_case '_VTF_data_corrected.mat'];
VTF_data_name=[defect_case '_VTF_data.mat'];

try
    load(VTF_data_name) % make sure same as in constraints function!!!
%     A=find(Pseudo_ductile_strain>1.6); % testing only
catch
    initial_data=[];initial_objective=[];initial_constraints=[];initial_errors=[];size_of_data=0;minInitS=[];minPdStrain=[];minUltStrain=[];minUltStrength=[];minYStrength=[];
    disp(['no ' defect_case ' VTF data - time to start from scratch!'])
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
     case {'eConstraint','eConstraintHybrid'}
%         initial_objective=-1*(wInitS*(Initial_stiffness)+wUltStrain*(Ultimate_strain)+wPdStrain*(Pseudo_ductile_strain)+wUltStrength*(Ultimate_strength)+wYStrength*(B));
        load('SOBO_FibreType_MinObesrved.mat')
        initial_objective=-1.*(wInitS.*(Initial_stiffness./minInitS)+wUltStrain.*(Ultimate_strain./minUltStrain)+wPdStrain.*(Pseudo_ductile_strain./minPdStrain)+wUltStrength.*(Ultimate_strength./minUltStrength)+wYStrength.*(B./minYStrength));
        minInitS=[];minPdStrain=[];minUltStrain=[];minUltStrength=[];minYStrength=[];
end

% set constraints and error traces
% initial_constraints=-50+Max_load_drop;
initial_constraints=define_initial_constraints_hybrid(VTF_data_name,e_constraint_name,e_constraint_min,n_e_constraint_slices,e_constraint_max,batch_index,defect_case);
initial_errors=2*isnan(initial_objective)-1;

size_of_data=size(Vc,1); % record size of dataset

% testing only!!!
% initial_data=initial_data(A,:);initial_objective=initial_objective(A);initial_constraints=initial_constraints(A,:);initial_errors=initial_errors(A);size_of_data=length(A);

end

