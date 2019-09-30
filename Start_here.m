%% A Bayesian optimisation routine for maximising the mechanical performance of aligned discontiuous composites
% current setup has a large amount of observed data, therefore will take
% some time (and RAM) to run. Minimum 16GB RAM recommended for use.
% James M. Finley, 10th September 2019, Imperial College London
function results = Start_here()

%% select primary objective objective function to optimise (select one from below)
% objective_function_1='Initial_stiffness';
objective_function_1='Ultimate_strain';
% objective_function_1='Pseudo_ductile_strain';
% objective_function_1='Ultimate_strength';
% objective_function_1='Yield_strength';


%% select secondary objective function to optimise (select one from below) - multi-objective optimisation is RAM intensive, therefore suggested to only use on a HPC
objective_function_2=''; % select this option single-objective optimisation of primary objective function
% objective_function_2='Initial_stiffness';
% objective_function_2='Ultimate_strain';
% objective_function_2='Pseudo_ductile_strain';
% objective_function_2='Ultimate_strength';
% objective_function_2='Yield_strength';


%% set optimisation parameters
number_of_iterations=1; % number of iterations of optimiser (or number of iterations per multi-objective potimisation slice if using multi-objective optimisation)
exploration_ratio=50; % exploration ratio (in percent) 

number_of_MOBO_slices=40; % number of multi-objective optimisation slices (different slices defined by the modified epsilon constraint method) - multi-objective optimisation only


%% configure the optimiser - DON'T MODIFY PAST THIS POINT
% check correct objective function names used
switch objective_function_1
    case{'Initial_stiffness'}
        objective_function_1_flag='wInitS';
    case{'Ultimate_strain'}
        objective_function_1_flag='wUltStrain';
    case{'Pseudo_ductile_strain'}
        objective_function_1_flag='wPdStrain';
    case{'Ultimate_strength'}
        objective_function_1_flag='wUltStrength';
    case{'Yield_strength'}
        objective_function_1_flag='wYStrength';
    otherwise
        error("incorrect objective function 1 name. Please choose from 'Initial_stiffness', 'Ultimate_strain', 'Pseudo_ductile_strain', 'Ultimate_strength', or 'Yield_strength'.");
end

single_multi_flag='MOBO'; % treat as multi-objective optimisation campaign unless otherwise stated
switch objective_function_2
    case ''
        single_multi_flag='SOBO'; % if no second objective function, treat as single-objective Bayesian optimisation (SOBO)
    case{'Initial_stiffness'} % set upper and lower constraint range
        e_constraint_min=0;
        e_constraint_max=6E+06;
    case {'Ultimate_strain','Pseudo_ductile_strain'}
        e_constraint_min=0;
        e_constraint_max=2000;
    case {'Ultimate_strength','Yield_strength'}
        e_constraint_min=0;
        e_constraint_max=4.5;
    otherwise
        error("incorrect objective function 2 name. Please choose from 'Initial_stiffness', 'Ultimate_strain', 'Pseudo_ductile_strain', 'Ultimate_strength', or 'Yield_strength', or for single-objective optimisation."); % error if wrong objective function name
end

%% run the optimiser 
switch single_multi_flag
    case 'SOBO' 
        results = Create_VTF_BayesianOptimizationObject(single_multi_flag,'single',objective_function_1,number_of_iterations,exploration_ratio,0,0,0,0,0,0,1); % run single-objective optimisation and store results
        gather_VTF_data('pristine') % update observed data
    case 'MOBO' %- very RAM intensive
        for ii=1:number_of_MOBO_slices % run a number of optimisations for multi-objective design case
            Create_VTF_BayesianOptimizationObject_eConstraint_hybrid(objective_function_2,e_constraint_min,number_of_MOBO_slices,e_constraint_max,objective_function_1_flag,single_multi_flag,'eConstraintHybrid',batch_index,number_of_iterations,exploration_ratio,0,0,0,0,0,0,1);
        end  
        % find results saved in your directory. 
        gather_VTF_data('pristine') % update observed data
        results=load('pristine_VTF_data.mat'); % load observed data
        figure()
        plot(results.(objective_function_2),results.(objective_function_1),'rx') % show observed data for the two selected objective functions
end

end