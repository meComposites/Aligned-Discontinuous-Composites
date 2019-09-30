function run_BayesianOptimisation_eConstraintHybrid(batch_index,defect_case,num_obj_function_evaluations,num_slices,exploration_ratio)

% set amount of defects for optimisation
switch defect_case
    case 'pristine'
        [pc,pInterface,pcFailed,pg,pgFailed,misalignment_switch]=deal(0);
    otherwise
        [pc,pInterface,pcFailed,pg,pgFailed]=deal(10);
        misalignment_switch=0;
end

% set up constrait names, maximise names, and constraint limits
constraint_names={'Yield_strength','Ultimate_strength','Initial_stiffness','Pseudo_ductile_strain','Ultimate_strain','Ultimate_strength'};
maximise_names={'wPdStrain','wUltStrain','wUltStrength','wYStrength','wUltStrength','wInitS'};
constraint_limits=[0 2000;0 2000;0 6E+5;0 4.5;0 4.5;0 2000];

% choose which parameter to optimise
pareto_index=ceil(batch_index/num_slices);

%set the slice index for the optimiser (which constraint limits to choose)
slice_index=mod(batch_index,num_slices);
if slice_index == 0
    slice_index=num_slices;
end

% run the optimisation
Create_VTF_BayesianOptimizationObject_eConstraint_hybrid(constraint_names{pareto_index},constraint_limits(pareto_index,1),num_slices,constraint_limits(pareto_index,2),maximise_names{pareto_index},'MOBO','eConstraintHybrid',slice_index,num_obj_function_evaluations,exploration_ratio,pc,pInterface,pcFailed,pg,pgFailed,misalignment_switch,1)

% write dummy file to workspace to tell HPC when it's successfully finished
fileID = fopen('sentinel.txt','w');
nbytes = fprintf(fileID,'optimisation complete');
end