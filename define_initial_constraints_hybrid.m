function initial_constraints=define_initial_constraints_hybrid(VTF_data_name,e_constraint_name,e_constraint_min,n_e_constraint_slices,e_constraint_max,batch_index,defect_case)

try
    load(VTF_data_name)
catch
    initial_data=[];initial_objective=[];initial_constraints=[];initial_errors=[];size_of_data=0;minInitS=[];minPdStrain=[];minUltStrain=[];minUltStrength=[];minYStrength=[];
    disp(['no ' defect_case ' VTF data - time to start from scratch!'])
    return
end

if batch_index == 1
    A=[];
else
    delta_constraint=(e_constraint_max-e_constraint_min)/n_e_constraint_slices;
    min_constraint=e_constraint_min+(batch_index-1)*delta_constraint;
    % max_constraint=e_constraint_min+(batch_index)*delta_constraint;
    
    % add little bit at end so that ultimate_strain case can get edges (doesn't
    % matter otherwise)
    A=eval(['min_constraint-' e_constraint_name '-0.001']);
    % B=eval([e_constraint_name '-max_constraint' '-0.001']);
end



% initial_constraints=[-50+Max_load_drop,A,B];
initial_constraints=[-50+Max_load_drop,A];
end