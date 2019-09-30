function initial_constraints=define_initial_constraints_noRBound(e_constraint_name,e_constraint_min,n_e_constraint_slices,batch_index,defect_case)

try
    load([defect_case '_VTF_data.mat'])
catch
    initial_data=[];initial_objective=[];initial_constraints=[];initial_errors=[];size_of_data=0;minInitS=[];minPdStrain=[];minUltStrain=[];minUltStrength=[];minYStrength=[];
    disp(['no ' defect_case ' VTF data - time to start from scratch!'])
    return
end

delta_constraint=(e_constraint_max-e_constraint_min)/n_e_constraint_slices;
min_constraint=e_constraint_min+(batch_index-1)*delta_constraint;

% add little bit at end so that ultimate_strain case can get edges (doesn't
% matter otherwise)
A=eval(['min_constraint-' e_constraint_name '-0.001']);
B=eval([e_constraint_name '-max_constraint' '-0.001']);

initial_constraints=[-50+Max_load_drop,A,B];
end