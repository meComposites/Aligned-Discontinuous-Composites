function [Objective_function_output,coupled_constraints] = VTF_objective_function_eConstraint_hybrid(single_multi_flag,objective_fun_flag,objective_name,e_constraint_name,e_constraint_min,n_e_constraint_slices,e_constraint_max,batch_index,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,wInitS,wUltStrain,wPdStrain,wUltStrength,wYStrength,pc,pInterface,pcFailed,pg,pgFailed,misalignment_switch,rng_switch,minInitS,minPdStrain,minUltStrain,minUltStrength,minYStrength)
% Objective function wrapper for VTF optimisation scheme
%   Takes input variables and evaluates the VTF at these points. Future
%   work will also save each data point to add to training database (for
%   multi-objective optimisation).

try 
    Specimen_dataset = Specimen_Microscale_FibreTypeBO(single_multi_flag,[objective_fun_flag '_max_' objective_name '_cons_' e_constraint_name '_' num2str(batch_index)],0.55,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,0.15,10,50000,80,80,pc,pInterface,pcFailed,pg,pInterface,pgFailed,misalignment_switch,rng_switch,0,1);
    
if isempty(Specimen_dataset.Yield_strength)==1
    Yield_strength=Specimen_dataset.Ultimate_strength
elseif isnan(Specimen_dataset.Yield_strength)==1
    Yield_strength=Specimen_dataset.Ultimate_strength
else
    Yield_strength=Specimen_dataset.Yield_strength
end

switch objective_fun_flag
    case {'eConstraint','eConstraintHybrid'}
        Objective_function_output=-1*(wInitS*(Specimen_dataset.Initial_stiffness)+wUltStrain*(Specimen_dataset.Ultimate_strain)+wPdStrain*(Specimen_dataset.Pseudo_ductile_strain)+wUltStrength*(Specimen_dataset.Ultimate_strength)+wYStrength*(Yield_strength));
    otherwise
        error('incorrect objective function flag');
end

% coupled_constraints=-50+Specimen_dataset.Max_load_drop;
Initial_stiffness=Specimen_dataset.Initial_stiffness
Ultimate_strain=Specimen_dataset.Ultimate_strain
Pseudo_ductile_strain=Specimen_dataset.Pseudo_ductile_strain
Ultimate_strength=Specimen_dataset.Ultimate_strength
e_constraint_value=eval(e_constraint_name); % now covers the case wher Yield_strength is NaN
coupled_constraints=determine_constraints_hybrid(e_constraint_value,e_constraint_min,n_e_constraint_slices,e_constraint_max,batch_index,Specimen_dataset.Max_load_drop);
load_drop_constraint=-50+Specimen_dataset.Max_load_drop;
if batch_index == 1
    min_constraint=[];
else
    min_constraint=coupled_constraints(:,2);
end
% max_constraint=coupled_constraints(:,3);

% coupled_constraints=[load_drop_constraint,min_constraint,max_constraint];
coupled_constraints=[load_drop_constraint,min_constraint]


catch e 
    disp(['error in file: ' e.stack(1).name ', line: ' num2str(e.stack(1).line) ': ']); 
    disp(e.message);
    
   Objective_function_output=NaN;
    % check whether this is okay
    coupled_constraints=NaN;
end

end

function coupled_constraints=determine_constraints_hybrid(e_constraint_value,e_constraint_min,n_e_constraint_slices,e_constraint_max,batch_index,Max_load_drop)

delta_constraint=(e_constraint_max-e_constraint_min)/n_e_constraint_slices;
min_constraint=e_constraint_min+(batch_index-1)*delta_constraint;
% max_constraint=e_constraint_min+(batch_index)*delta_constraint;

% add little bit at end so that ultimate_strain case can get edges (doesn't
% matter otherwise)
A=min_constraint- e_constraint_value -0.001;
% B=e_constraint_value -max_constraint-0.001;

% coupled_constraints=[-50+Max_load_drop,A,B];
coupled_constraints=[-50+Max_load_drop,A];
end
