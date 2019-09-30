function investigate_MO_pareto_variability(batch_index,defect_case)

load([defect_case '_MO_pareto_variability_inputs.mat']);

run_setting=mod(batch_index,num_files);
if run_setting==0
    run_setting=num_files;
end

repeat_counter=ceil(batch_index/num_files);

thisFilename=char(Filename(run_setting));
thisFilename=thisFilename(1:strfind(thisFilename,'.mat')-1);

Specimen_dataset = Specimen_Microscale_FibreTypeBO('MO_variability',[char(MO_case(run_setting)) '_' thisFilename '_' num2str(repeat_counter)],0.55,Carbon_fibre(run_setting),Glass_fibre(run_setting),Vc(run_setting),lf(run_setting)*1000,SmAvg(run_setting),G(run_setting),GiicmAvg(run_setting),0.15,10,50000,80,80,pc_fragmented(run_setting),p_interface(run_setting),pc_vacancy(run_setting),pg_fragmented(run_setting),p_interface(run_setting),pg_vacancy(run_setting),misalignment_case(run_setting),1,0,1);

% write dummy file to workspace to tell HPC when it's successfully finished
fileID = fopen('sentinel.txt','w');
nbytes = fprintf(fileID,'optimisation complete');

end