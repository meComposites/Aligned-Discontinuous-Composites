%% creates estimated mean (and estimated stdDev) outputs from the surrogate model, for the input data supplied. Must use up-to-date surrogate model, fitted to all available VTF data.
...Similar to sensitivity_analysis, but only for estimated mean output, not for modifications to the intended inputs
    
% Need to add SOBO results to path

function estimate_VTF_data(defect_case,data_trace_case)
dbstop if error
propertyIndices=[1,3,4,5];
% filenameIndices={'max_wInitS_contrained_Ultimate_strength';'max_wUltStrain_contrained_Ultimate_strength';'max_wPdStrain_contrained_Yield_strength';'max_wYStrength_contrained_Pseudo_ductile_strain';'max_wUltStrength_contrained_Initial_stiffness'};
filenameIndices={'max_wInitS_contrained_Ultimate_strength';'max_wUltStrain_contrained_Ultimate_strength';'max_wPdStrain_contrained_Yield_strength';'max_wYStrength_contrained_Pseudo_ductile_strain';'max_wUltStrength_contrained_Initial_stiffness'};
load(['SOBO_' defect_case '_FibreType_MinObesrved.mat']);
minObservedValues=[minInitS,minUltStrain,minPdStrain,minUltStrength,minYStrength];
switch defect_case
    case 'pristine'
        switch data_trace_case
            case 'VTF_data_trace' % normal VTF data for estimated mean etc.
%                 load(['results-second_pass/corrected_models/MOBO_eConstraintHybrid_' defect_case '_results_' filenameIndices{1} '_segment_1.mat'])
                load(['results-second_pass/corrected_models/non-corrected_normalised_models/MOBO_eConstraintHybrid_' defect_case '_results_' filenameIndices{1} '_segment_1.mat'])
                XTrace=results.XTrace;
            case 'corrected_VTF_data_trace'
                %         load(['pristine_VTF_data_corrected_sensitivityInputs.mat'])
                %         lf=lf*1000;
                %         XTrace=table(Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg);
            otherwise
                % %         load(['D:/Box Sync/HiPerDuCT PhD/Modelling/Optimisation/results/MOBO/additive/pristine/plots/paretoFronts/' defect_case '_' data_trace_case '_pareto_inputs.mat']) % make sure data_trace_case contains filename you want!
                % %         load(['matlab_scripts/' defect_case '_observedParetoFront_inputs_alternative.mat']);
                %         load(['matlab_scripts/' defect_case '_observedParetoFront_inputs_alternative_corrected.mat']);
                %         XTrace=table(Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg); % check lf is consistent
                %         MO_case=MO_case;
        end
    case 'defective'
        switch data_trace_case
            case 'VTF_data_trace' % normal VTF data for estimated mean etc.
                load(['SOBO_single_' defect_case '_results_' filenameIndices{1} '.mat'])
                XTrace=results.XTrace;
        end
end

for propertyIndex=1:5
    switch defect_case
        case 'pristine'
            load(['results-second_pass/corrected_models/non-corrected_normalised_models/MOBO_eConstraintHybrid_' defect_case '_results_' filenameIndices{propertyIndices(propertyIndex)} '_segment_1.mat'])
            
        case 'defective'
            load(['SOBO_single_' defect_case '_results_' filenameIndices{propertyIndex} '.mat'])
    end
             
    %     testin only
    %     load(['pristine_VTF_data_corrected_sensitivityInputs.mat'])
    %     XTrace=[XTrace;table(Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg)];
    
    if propertyIndex==1
        [meanProperties,stdDevProperties,constraints,stdDevConstraints,constraintViolations]=initialise_data(size(XTrace,1));
    end
    
    
    
    
    [meanProperties(:,propertyIndex),stdDevProperties(:,propertyIndex)] = predictObjective(results,XTrace);
    [constraints(:,propertyIndex),stdDevConstraints(:,propertyIndex)] = predictConstraints(results,XTrace);
    
    switch defect_case
        case 'pristine'
            meanProperties(:,propertyIndex)=meanProperties(:,propertyIndex).*-1;
%             stdDevProperties(:,propertyIndex)=stdDevProperties(:,propertyIndex).*minObservedValues(propertyIndex);
        case 'defective'
            meanProperties(:,propertyIndex)=meanProperties(:,propertyIndex).*-1;
    end
    
end

% StdProperties=meanProperties-1.96.*stdDevProperties; % based on normal distribution with 95% Std
stdProperties=meanProperties-1.*stdDevProperties;

constraintViolations=reshape(mode(constraints<=0,2).*(-2)+1,size(XTrace,1),1); % constraint violations taken as moe for all 5 properties

switch data_trace_case
    case 'VTF_data_trace'
        saveName=['matlab_scripts-second_pass/' defect_case '_estimated_VTF_data.mat'];
        save(saveName,'meanProperties','stdDevProperties','stdProperties','constraints','stdDevConstraints','constraintViolations','-v7.3');
        save(['results-second_pass/estimatedParetoFronts/' defect_case '_estimated_VTF_data.mat'],'meanProperties','stdDevProperties','stdProperties','constraints','stdDevConstraints','constraintViolations','-v7.3');
    otherwise
        saveName=['matlab_scripts-second_pass/' defect_case '_' data_trace_case '_variability_analysis.mat'];
        save(saveName,'meanProperties','stdDevProperties','stdProperties','constraints','stdDevConstraints','constraintViolations','-v7.3');
        save(['results-second_pass/estimatedParetoFronts/' defect_case '_' data_trace_case '_variability_analysis.mat'],'meanProperties','stdDevProperties','stdProperties','constraints','stdDevConstraints','constraintViolations','-v7.3');
end
%
% save(saveName,'meanProperties','stdDevProperties','stdProperties','constraints','stdDevConstraints','constraintViolations',...
%     'meanMeanProperties','stdDevMeanProperties','stdMeanProperties','meanStdProperties','stdDevStdProperties','stdStdProperties',...
%     'CIMeanProperties','minMeanProperties','minMeanSensitivitySortIndex','minStdProperties','minStdSensitivitySortIndex','-v7.3');

end

function [meanProperties,stdDevProperties,constraints,stdDevConstraints,constraintViolations]=initialise_data(resultsSize)
[meanProperties,stdDevProperties,constraints,stdDevConstraints]=deal(nan(resultsSize,5));
constraintViolations=nan(resultsSize,1);
end
