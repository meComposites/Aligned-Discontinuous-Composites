function [paretoOptimal_ultStrength_initS,paretoOptimal_ultStrength_ultStrain,paretoOptimal_yieldStrength_pdStrain,...
    paretoOptimal_ultStrength_initS_filenames,paretoOptimal_ultStrength_ultStrain_filenames,paretoOptimal_yieldStrength_pdStrain_filenames,...
    paretoOptimal_ultStrength_initS_fileIndices,paretoOptimal_ultStrength_ultStrain_fileIndices,paretoOptimal_yieldStrength_pdStrain_fileIndices]=...
    find_pareto_optimals(defect_case,observed_estimated_case,best_worse_case,Initial_data,varargin)

dbstop if error
maxPermittedLoadDrop=50;
MO_case_check={'ultStrength_initS','ultStrength_ultStrain','yieldStrength_pdStrain'};

if nargin == 4
    Filename=strings(size(Initial_dataset,1),1);
else
    % load(['matlab_scripts/' defect_case '_VTF_data_corrected.mat']);
    % load data, i.e. performance metrics, filenames, and constraint
    % violations
    switch observed_estimated_case
        case 'observed'
%             load(['matlab_scripts/' defect_case '_VTF_data.mat']);
            load(['matlab_scripts-second_pass/' defect_case '_VTF_data_corrected.mat']);
            % properties=["initialStiffness","ultimateStrain","pDStrain","ultimateStrength","yieldStrength"];
            Initial_dataset=[Initial_stiffness,Ultimate_strain,Pseudo_ductile_strain,Ultimate_strength,Yield_strength];
            constraintViolations=(Max_load_drop<maxPermittedLoadDrop)*-2+1;
        case 'estimatedMean'
            switch defect_case
                case 'defective'
                    load(['matlab_scripts/' defect_case '_estimated_VTF_data.mat']);
                otherwise
                    load(['results-second_pass/estimatedParetoFronts/' defect_case '_estimated_VTF_data_corrected.mat']);
            end
            Initial_dataset=[meanProperties(:,1),meanProperties(:,2),meanProperties(:,3),meanProperties(:,4),meanProperties(:,5)];
            constraintViolations=mode(constraints<0,2)*-2+1;
            Filename=strings(size(meanProperties,1),1);
        case 'initial_sensitivity_Pareto'
            load('results-second_pass/sensitivityAnalysis/initial_sensitivity_analysis/paretoFrontSensitivityAnalysis.mat')
            Initial_dataset=[meanProperties(:,:,1);meanProperties(:,:,2);meanProperties(:,:,3);meanProperties(:,:,4);meanProperties(:,:,5);meanProperties(:,:,7);meanProperties(:,:,8);meanProperties(:,:,9);meanProperties(:,:,10);meanProperties(:,:,11)];
            constraints=[constraints(:,:,1);constraints(:,:,2);constraints(:,:,3);constraints(:,:,4);constraints(:,:,5);constraints(:,:,7);constraints(:,:,8);constraints(:,:,9);constraints(:,:,10);constraints(:,:,11)];
            constraintViolations=mode(constraints<0,2)*-2+1;
            Filename=[sensitivityCase(:,1);sensitivityCase(:,2);sensitivityCase(:,3);sensitivityCase(:,4);sensitivityCase(:,5);sensitivityCase(:,7);sensitivityCase(:,8);sensitivityCase(:,9);sensitivityCase(:,10);sensitivityCase(:,11)];
            MO_case=repmat(MO_case,10,1);
        case 'no_small_hybrids_sensitivity_Pareto'
            load('results-second_pass/sensitivityAnalysis/no_small_hybrids/paretoFrontSensitivityAnalysis.mat')
            Initial_dataset=[meanProperties(:,:,1);meanProperties(:,:,2);meanProperties(:,:,3);meanProperties(:,:,4);meanProperties(:,:,5);meanProperties(:,:,7);meanProperties(:,:,8);meanProperties(:,:,9);meanProperties(:,:,10);meanProperties(:,:,11)];
            constraints=[constraints(:,:,1);constraints(:,:,2);constraints(:,:,3);constraints(:,:,4);constraints(:,:,5);constraints(:,:,7);constraints(:,:,8);constraints(:,:,9);constraints(:,:,10);constraints(:,:,11)];
            constraintViolations=mode(constraints<0,2)*-2+1;
            Filename=[sensitivityCase(:,1);sensitivityCase(:,2);sensitivityCase(:,3);sensitivityCase(:,4);sensitivityCase(:,5);sensitivityCase(:,7);sensitivityCase(:,8);sensitivityCase(:,9);sensitivityCase(:,10);sensitivityCase(:,11)];
            MO_case=repmat(MO_case,10,1);
%         case 'estimatedMean_fromSensitivityAnalysis'
%             load(['results-second_pass/estimatedParetoFronts/' defect_case '_estimated_VTF_data_corrected.mat']);
%             load('results-second_pass/sensitivityAnalysis/initial_sensitivity_analysis/local_10pc_pristine_corrected_VTF_data_trace_sensitivity_analysis.mat')
%             Initial_dataset=[meanProperties(:,1,6),meanProperties(:,2,6),meanProperties(:,3,6),meanProperties(:,4,6),meanProperties(:,5,6)];
%             constraintViolations=mode(reshape(constraints(:,:,6),size(constraints,1),size(constraints,2))<0,2)*-2+1;
%             Filename=strings(size(meanProperties,1),1);
%         case 'estimatedStd'
%             load(['matlab_scripts/' defect_case '_estimated_pareto_data.mat']);
%             Initial_dataset=[stdProperties(:,1),stdProperties(:,2),stdProperties(:,3),stdProperties(:,4),stdProperties(:,5)];
%             constraintViolations=mode(constraints<0,2)*-2+1;
%             Filename=strings(size(stdProperties,1),1);
%             
%             % for sensitivity cases, remove the 6th entry, as it is the
%             % baseline
%         case 'sensitivityMean'
%             %             load(['matlab_scripts/local_10pc_' defect_case '_sensitivity_analysis.mat']);
%             sensitivityMatFile=matfile(['matlab_scripts/local_10pc_' defect_case '_sensitivity_analysis.mat']);
%             meanProperties=sensitivityMatFile.meanProperties;
%             constraints=sensitivityMatFile.constraints;
%             constraints=[constraints(:,:,1);constraints(:,:,2);constraints(:,:,3);constraints(:,:,4);constraints(:,:,5);constraints(:,:,7);constraints(:,:,8);constraints(:,:,9);constraints(:,:,10);constraints(:,:,11)];
%             constraintViolations=mode(constraints<0,2)*-2+1;
%             sensitivityCase=sensitivityMatFile.sensitivityCase;
%             Initial_dataset=[meanProperties(:,:,1);meanProperties(:,:,2);meanProperties(:,:,3);meanProperties(:,:,4);meanProperties(:,:,5);meanProperties(:,:,7);meanProperties(:,:,8);meanProperties(:,:,9);meanProperties(:,:,10);meanProperties(:,:,11)];
%             Filename=[sensitivityCase(:,1);sensitivityCase(:,2);sensitivityCase(:,3);sensitivityCase(:,4);sensitivityCase(:,5);sensitivityCase(:,7);sensitivityCase(:,8);sensitivityCase(:,9);sensitivityCase(:,10);sensitivityCase(:,11)];
%         case 'sensitivityStd'
%             %             load(['matlab_scripts/local_10pc_' defect_case '_sensitivity_analysis.mat']);
%             sensitivityMatFile=matfile(['matlab_scripts/local_10pc_' defect_case '_sensitivity_analysis.mat']);
%             stdProperties=sensitivityMatFile.stdProperties;
%             constraints=sensitivityMatFile.constraints;
%             constraints=[constraints(:,:,1);constraints(:,:,2);constraints(:,:,3);constraints(:,:,4);constraints(:,:,5);constraints(:,:,7);constraints(:,:,8);constraints(:,:,9);constraints(:,:,10);constraints(:,:,11)];
%             constraintViolations=mode(constraints<0,2)*-2+1;
%             sensitivityCase=sensitivityMatFile.sensitivityCase;
%             Initial_dataset=[stdProperties(:,:,1);stdProperties(:,:,2);stdProperties(:,:,3);stdProperties(:,:,4);stdProperties(:,:,5);stdProperties(:,:,7);stdProperties(:,:,8);stdProperties(:,:,9);stdProperties(:,:,10);stdProperties(:,:,11)];
%             Filename=[sensitivityCase(:,1);sensitivityCase(:,2);sensitivityCase(:,3);sensitivityCase(:,4);sensitivityCase(:,5);sensitivityCase(:,7);sensitivityCase(:,8);sensitivityCase(:,9);sensitivityCase(:,10);sensitivityCase(:,11)];
%         case 'paretoSensitivityMean'
% %             sensitivityMatFile=matfile(['matlab_scripts/' defect_case '_observedParetoFront_inputs_alternative_sensitivity_analysis.mat']);
% %             sensitivityMatFile=matfile('local_10pc_pristine_observedParetoFronts_corrected_modsSameFibre_sensitivity_analysis.mat');
%             sensitivityMatFile=matfile('alternative_paretoFrontSensitivityAnalysis_worst.mat');
%             meanProperties=sensitivityMatFile.meanProperties;
%             constraints=sensitivityMatFile.constraints;
%             constraints=[constraints(:,:,1);constraints(:,:,2);constraints(:,:,3);constraints(:,:,4);constraints(:,:,5);constraints(:,:,7);constraints(:,:,8);constraints(:,:,9);constraints(:,:,10);constraints(:,:,11)];
%             constraintViolations=mode(constraints<0,2)*-2+1;
%             sensitivityCase=sensitivityMatFile.sensitivityCase;
%             Initial_dataset=[meanProperties(:,:,1);meanProperties(:,:,2);meanProperties(:,:,3);meanProperties(:,:,4);meanProperties(:,:,5);meanProperties(:,:,7);meanProperties(:,:,8);meanProperties(:,:,9);meanProperties(:,:,10);meanProperties(:,:,11)];
%             Filename=[sensitivityCase(:,1);sensitivityCase(:,2);sensitivityCase(:,3);sensitivityCase(:,4);sensitivityCase(:,5);sensitivityCase(:,7);sensitivityCase(:,8);sensitivityCase(:,9);sensitivityCase(:,10);sensitivityCase(:,11)];
%             MO_case=repmat(sensitivityMatFile.MO_case,10,1);
%         case 'paretoSensitivityStd'
%             sensitivityMatFile=matfile(['matlab_scripts/local_10pc_' defect_case '_observedParetoFronts_sensitivity_analysis.mat']);
%             stdProperties=sensitivityMatFile.stdProperties;
%             constraints=sensitivityMatFile.constraints;
%             constraints=[constraints(:,:,1);constraints(:,:,2);constraints(:,:,3);constraints(:,:,4);constraints(:,:,5);constraints(:,:,7);constraints(:,:,8);constraints(:,:,9);constraints(:,:,10);constraints(:,:,11)];
%             constraintViolations=mode(constraints<0,2)*-2+1;
%             sensitivityCase=sensitivityMatFile.sensitivityCase;
%             Initial_dataset=[stdProperties(:,:,1);stdProperties(:,:,2);stdProperties(:,:,3);stdProperties(:,:,4);stdProperties(:,:,5);stdProperties(:,:,7);stdProperties(:,:,8);stdProperties(:,:,9);stdProperties(:,:,10);stdProperties(:,:,11)];
%             Filename=[sensitivityCase(:,1);sensitivityCase(:,2);sensitivityCase(:,3);sensitivityCase(:,4);sensitivityCase(:,5);sensitivityCase(:,7);sensitivityCase(:,8);sensitivityCase(:,9);sensitivityCase(:,10);sensitivityCase(:,11)];
%             MO_case=repmat(sensitivityMatFile.MO_case,10,1);
%         case 'estimatedMean_correctedTest'
%             load(['alternative_paretoFrontSensitivityAnalysis_worst_correctedTest.mat']);
%             meanProperties=reshape(meanProperties(:,:,6),size(meanProperties,1),size(meanProperties,2));
%             constraints=reshape(constraints(:,:,6),size(meanProperties,1),size(meanProperties,2));
%             Initial_dataset=[meanProperties(:,1),meanProperties(:,2),meanProperties(:,3),meanProperties(:,4),meanProperties(:,5)];
%             constraintViolations=mode(constraints<0,2)*-2+1;
%             Filename=strings(size(meanProperties,1),1);
%         case 'paretoSensitivityMean_correctedTest'
% %             sensitivityMatFile=matfile(['matlab_scripts/' defect_case '_observedParetoFront_inputs_alternative_sensitivity_analysis.mat']);
% %             sensitivityMatFile=matfile('local_10pc_pristine_observedParetoFronts_corrected_modsSameFibre_sensitivity_analysis.mat');
%             sensitivityMatFile=matfile('alternative_paretoFrontSensitivityAnalysis_worst_correctedTest.mat');
%             meanProperties=sensitivityMatFile.meanProperties;
%             constraints=sensitivityMatFile.constraints;
%             constraints=[constraints(:,:,1);constraints(:,:,2);constraints(:,:,3);constraints(:,:,4);constraints(:,:,5);constraints(:,:,7);constraints(:,:,8);constraints(:,:,9);constraints(:,:,10);constraints(:,:,11)];
%             constraintViolations=mode(constraints<0,2)*-2+1;
%             sensitivityCase=sensitivityMatFile.sensitivityCase;
%             Initial_dataset=[meanProperties(:,:,1);meanProperties(:,:,2);meanProperties(:,:,3);meanProperties(:,:,4);meanProperties(:,:,5);meanProperties(:,:,7);meanProperties(:,:,8);meanProperties(:,:,9);meanProperties(:,:,10);meanProperties(:,:,11)];
%             Filename=[sensitivityCase(:,1);sensitivityCase(:,2);sensitivityCase(:,3);sensitivityCase(:,4);sensitivityCase(:,5);sensitivityCase(:,7);sensitivityCase(:,8);sensitivityCase(:,9);sensitivityCase(:,10);sensitivityCase(:,11)];
%             MO_case=repmat(sensitivityMatFile.MO_case,10,1);
%         case 'estimatedMean_correctedTest_noSmallHybrids'
%             load(['local_10pc_pristine_sensitivity_analysis_correctedTest_noSmallHybrids.mat']);
%             meanProperties=reshape(meanProperties(:,:,6),size(meanProperties,1),size(meanProperties,2));
%             constraints=reshape(constraints(:,:,6),size(meanProperties,1),size(meanProperties,2));
%             Initial_dataset=[meanProperties(:,1),meanProperties(:,2),meanProperties(:,3),meanProperties(:,4),meanProperties(:,5)];
%             constraintViolations=mode(constraints<0,2)*-2+1;
%             Filename=strings(size(meanProperties,1),1);
%         case 'paretoSensitivityMean_correctedTest_noSmallHybrids'
%             load(['alternative_paretoFrontSensitivityAnalysis_worst_correctedTest_noSmallHybrids.mat']);
%             sensitivityMatFile=matfile('alternative_paretoFrontSensitivityAnalysis_worst_correctedTest_noSmallHybrids.mat');
%             meanProperties=sensitivityMatFile.meanProperties;
%             constraints=sensitivityMatFile.constraints;
%             constraints=[constraints(:,:,1);constraints(:,:,2);constraints(:,:,3);constraints(:,:,4);constraints(:,:,5);constraints(:,:,7);constraints(:,:,8);constraints(:,:,9);constraints(:,:,10);constraints(:,:,11)];
%             constraintViolations=mode(constraints<0,2)*-2+1;
%             sensitivityCase=sensitivityMatFile.sensitivityCase;
%             Initial_dataset=[meanProperties(:,:,1);meanProperties(:,:,2);meanProperties(:,:,3);meanProperties(:,:,4);meanProperties(:,:,5);meanProperties(:,:,7);meanProperties(:,:,8);meanProperties(:,:,9);meanProperties(:,:,10);meanProperties(:,:,11)];
%             Filename=[sensitivityCase(:,1);sensitivityCase(:,2);sensitivityCase(:,3);sensitivityCase(:,4);sensitivityCase(:,5);sensitivityCase(:,7);sensitivityCase(:,8);sensitivityCase(:,9);sensitivityCase(:,10);sensitivityCase(:,11)];
%             MO_case=repmat(sensitivityMatFile.MO_case,10,1);
    end
    
end

% only use data where the load drop constraint is satisfied
fileIndices=[1:size(Initial_dataset,1)]';
Initial_dataset=Initial_dataset(constraintViolations<0,:);
Filename=Filename(constraintViolations<0);
fileIndices=fileIndices(constraintViolations<0);

switch observed_estimated_case
        case {'paretoSensitivityMean','paretoSensitivityStd','paretoSensitivityMean_correctedTest','paretoSensitivityMean_correctedTest_noSmallHybrids'}
           MO_case=MO_case(constraintViolations<0);
           Initial_dataset_copy=Initial_dataset;
           Filename_copy=Filename;
           fileIndices_copy=fileIndices;
end

% for loop to do max stiffness & strength, max ult strain and ult strength,
% and max yield strength and pdstrain
f_index=[4,4,5];
g_index=[1,2,3];

results_collection=cell(3,3);

for ii=1:3
    switch observed_estimated_case
        case {'paretoSensitivityMean','paretoSensitivityStd','paretoSensitivityMean_correctedTest','paretoSensitivityMean_correctedTest_noSmallHybrids'}
%             Initial_dataset=Initial_dataset_copy;
            Initial_dataset=Initial_dataset_copy(strcmp(MO_case,MO_case_check{ii})==1,:);
            Filename=Filename_copy(strcmp(MO_case,MO_case_check{ii})==1,:);
            fileIndices=fileIndices_copy(strcmp(MO_case,MO_case_check{ii})==1,:);
    end
    
    % find working datasets for each of the optimality metrics
    switch best_worse_case
        case 'worst'
            f_x=-1.*Initial_dataset(:,f_index(ii));
            g_x=-1.*Initial_dataset(:,g_index(ii));
        otherwise
            f_x=Initial_dataset(:,f_index(ii));
            g_x=Initial_dataset(:,g_index(ii));
    end
    
    working_Filename=cellstr(Filename);
    working_fileIndex=fileIndices;
    
    optimal_g=zeros(size(g_x,1),1);
    optimal_f=zeros(size(g_x,1),1);
    optimal_filename=cell(size(g_x,1),1);
    optimal_fileIndex=NaN(size(g_x,1),1);
    optimal_index=1;
    
    while size(g_x,1)>0
        % find the maximum value of g_x in the working dataset
        [g_x_max,g_x_max_index]=max(g_x);
        
        % if more than one max, then throw error (for now)
        if size(g_x_max,1)>1
            error('multiple maxima');
        end
        
        % add the maximum to the optimal set, and record the filename
        optimal_g(optimal_index)=g_x_max;
        optimal_f(optimal_index)=f_x(g_x_max_index);
        optimal_filename{optimal_index}=working_Filename{g_x_max_index};
        optimal_fileIndex(optimal_index)=working_fileIndex(g_x_max_index);
        
        % remove any sub-optimal members of the working set (f_x must be
        % larger than current optimal_f_x value)
        g_x=g_x(f_x>optimal_f(optimal_index));
        working_Filename=working_Filename(f_x>optimal_f(optimal_index));
        working_fileIndex=working_fileIndex(f_x>optimal_f(optimal_index));
        f_x=f_x(f_x>optimal_f(optimal_index)); % this one goes last to prevent errors?
        
        % increment the optimal dataset index for the next optimal point
        optimal_index=optimal_index+1;
    end
    
    % get rid of redundant zeros in optimal_f (depends on whether you are
    % looking at best or worst case
    switch best_worse_case
        case 'worst'
            optimal_f=optimal_f(optimal_f<0);
        otherwise
            optimal_f=optimal_f(optimal_f>0);
    end    
    optimal_g=optimal_g(1:size(optimal_f,1));
    optimal_filename=string(optimal_filename(1:size(optimal_g,1)));
    optimal_fileIndex=optimal_fileIndex(1:size(optimal_f,1));
    
    results_collection{1,ii}=optimal_g;
    results_collection{2,ii}=optimal_f;
    results_collection{3,ii}=optimal_filename;
    results_collection{4,ii}=optimal_fileIndex;
end

switch best_worse_case
    case 'worst'
        for ii=1:size(results_collection,1)-2
            for jj=1:size(results_collection,2)
                results_collection{ii,jj}=-1.*results_collection{ii,jj};
            end
        end
end

paretoOptimal_ultStrength_initS=[results_collection{1,1},results_collection{2,1}];
paretoOptimal_ultStrength_initS_filenames=results_collection{3,1};
paretoOptimal_ultStrength_initS_fileIndices=results_collection{4,1};

paretoOptimal_ultStrength_ultStrain=[results_collection{1,2},results_collection{2,2}];
paretoOptimal_ultStrength_ultStrain_filenames=results_collection{3,2};
paretoOptimal_ultStrength_ultStrain_fileIndices=results_collection{4,2};

paretoOptimal_yieldStrength_pdStrain=[results_collection{1,3},results_collection{2,3}];
paretoOptimal_yieldStrength_pdStrain_filenames=results_collection{3,3};
paretoOptimal_yieldStrength_pdStrain_fileIndices=results_collection{4,3};

end