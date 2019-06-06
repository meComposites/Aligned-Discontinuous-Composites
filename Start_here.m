%% Wrapper function for Virtual testing Framework, used for journal paper titled "The influence of variability and defects on the mechanical performance of tailorable composites". Code created by James M. Finley, Joel Henry, and Soraia Pimenta.
% Press play to run a working example of the code. Inputs can be modified
% here too (see below). More detailed modifications to the code (i.e. to
% modify fibre-type arrangements can be found in
% "Specimen_Microscale_EoR_Defects_FailureModes_General.m", while plotting
% commands for the het maps can be found in
% "PostProcess_EffectsOfRandomness_Results.m"


function Specimen_dataset = Start_here()
% running this code will create an object called specimen_dataset, which
% contains all data about the virtual specimen, including stress-strain
% curves etc. (see "Specimen.m" for details on the class definition and
% properties of the virtual specimen). A .mat file is also created
% conatining the specimen_dataset object (for running on a HPC cluster).

repeat_counter=1; % For when running repeat runs, to give unique ID to results

material_system=1; % 1 = HMC/EG/EP (intermingled), 2 = HSC/EP, 3 = HSC/PP

lf=3000; % fibre length (microns)
lspecimen=50000; % specimen length (microns)

n_rows=80; % no. fibres in y-dimension of cross-section
n_columns=80; % no. fibres in x-direction of cross-section

pc_fragmented = 0; % percetange of fragmented carbon fibres (0 - 100 %)
pc_vacancy = 0; % percentage of carbon fibres with fibre vacancy defects (0 - 100 %) 
p_interface = 0; %percentage of defective interactions (0 - 100 %)
pg_fragmented = 0; % percetange of fragmented glass fibres (0 - 100 %)
pg_vacancy = 0; % percentage of glass fibres with fibre vacancy defects (0 - 100 %) 

misalignment_case = 0; % 0 = no misalignment, 1 = HiPerDiF misalignment, 2 = Sanadi misalignment (see Figure 2b)

rng_switch = 1; % 0 = random number generator starts with repeat counter setting, 1 = mersenne twister rng (for HPC), otherwise, rng(other number)

failure_modes_switch = 0; % 1 = output details for damage events (used for damage event heat maps) 0 = suppress output

parallelisation_switch = 1; % 1 = use multiple CPUs for parallel calculation of RVEs 0 = serial calculation of RVEs

% Run the VTF
Specimen_dataset = Specimen_Microscale_EoR_Defects_FailureModes_General(repeat_counter,material_system,lf,lspecimen,n_rows,n_columns,pc_fragmented,p_interface,pc_vacancy,pg_fragmented,p_interface,pg_vacancy,misalignment_case,rng_switch,failure_modes_switch,parallelisation_switch);

end