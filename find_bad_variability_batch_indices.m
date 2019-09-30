function [bad_batch_indices,bad_batch_filenames]=find_bad_variability_batch_indices(batch_size,defect_case,save_flag)

% find file data in folder
fileData=dir(['G:\Optimisation\VTF_run_data-second_pass\variability\*.mat']);

% load input data
load([defect_case '_MO_pareto_variability_inputs.mat']);

bad_batch_indices=NaN(batch_size,1);
bad_batch_filenames=cell(batch_size,1);

% loop through for each batch run on the HPC
for batch_index=1:batch_size
    
    % find run setting an repeat counter
    run_setting=mod(batch_index,num_files);
    if run_setting==0
        run_setting=num_files;
    end
    
    repeat_counter=ceil(batch_index/num_files);
    
    thisFilename=char(Filename(run_setting));
    thisFilename=thisFilename(1:strfind(thisFilename,'.mat')-1);
    
    Filepath=['VTF_MO_variability_' char(MO_case(run_setting)) '_' thisFilename '_' num2str(repeat_counter) '_' defect_case '_run_.mat'];
    
    isFound=zeros(size(fileData,1),1);
    for ii=1:size(fileData,1)
        isFound(ii)=strcmp(Filepath,fileData(ii).name);
    end
    
    
    if sum((isFound))>0
        if sum((isFound))>1
            error('repeated data')
        else
            continue
        end
    else
        bad_batch_indices(batch_index)=batch_index;
        bad_batch_filenames{batch_index}=Filepath;
    end
    
end

bad_batch_indices=bad_batch_indices(isnan(bad_batch_indices)~=1);
bad_batch_filenames = bad_batch_filenames(~cellfun('isempty',bad_batch_filenames));

switch save_flag
    case 'save'
        save(['matlab_scripts-second_pass/' defect_case '_bad_variability_batch_data.mat'],'bad_batch_indices','bad_batch_filenames','-v7.3');
    otherwise
        
end
end