function [Objective_function_output,load_drop_constraint] = VTF_objective_function(single_multi_flag,objective_fun_flag,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,wInitS,wUltStrain,wPdStrain,wUltStrength,wYStrength,pc,pInterface,pcFailed,pg,pgFailed,misalignment_switch,rng_switch,minInitS,minPdStrain,minUltStrain,minUltStrength,minYStrength)
% Objective function wrapper for VTF optimisation scheme
%   Takes input variables and evaluates the VTF at these points. Future
%   work will also save each data point to add to training database (for
%   multi-objective optimisation).

try 
    Specimen_dataset = Specimen_Microscale_FibreTypeBO(single_multi_flag,objective_fun_flag,0.55,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,0.15,10,50000,80,80,pc,pInterface,pcFailed,pg,pInterface,pgFailed,misalignment_switch,rng_switch,0,0);
    
if isempty(Specimen_dataset.Yield_strength)==1
    Yield_strength=Specimen_dataset.Ultimate_strength;
elseif isnan(Specimen_dataset.Yield_strength)==1
    Yield_strength=Specimen_dataset.Ultimate_strength;
else
    Yield_strength=Specimen_dataset.Yield_strength;
end

switch objective_fun_flag
    case 'single'
        Objective_function_output=-1*(wInitS*(Specimen_dataset.Initial_stiffness)+wUltStrain*(Specimen_dataset.Ultimate_strain)+wPdStrain*(Specimen_dataset.Pseudo_ductile_strain)+wUltStrength*(Specimen_dataset.Ultimate_strength)+wYStrength*(Yield_strength));
    case 'multiplicative'
        Objective_function_output=-1*(wInitS*(Specimen_dataset.Initial_stiffness)*wUltStrength*(Specimen_dataset.Ultimate_strength)+wUltStrength*(Specimen_dataset.Ultimate_strength)*wUltStrain*(Specimen_dataset.Ultimate_strain)+wYStrength*(Yield_strength)*wPdStrain*(Specimen_dataset.Pseudo_ductile_strain));
    case {'additive','additive_detailed_ultStrengthUltStrain','additive_detailed_yStrengthPdStrain'}
        Objective_function_output=-1*(wInitS*(Specimen_dataset.Initial_stiffness/minInitS)+wUltStrain*(Specimen_dataset.Ultimate_strain/minUltStrain)+wPdStrain*(Specimen_dataset.Pseudo_ductile_strain/minPdStrain)+wUltStrength*(Specimen_dataset.Ultimate_strength/minUltStrength)+wYStrength*(Yield_strength/minYStrength));
end

load_drop_constraint=-50+Specimen_dataset.Max_load_drop;

catch e 
    disp(['error in file: ' e.stack(1).name ', line: ' num2str(e.stack(1).line) ': ']); 
    disp(e.message);
    
   Objective_function_output=NaN;
    % check whether this is okay
    load_drop_constraint=NaN;
end

end

