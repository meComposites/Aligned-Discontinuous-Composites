function sensitivity_analysis(defect_case,data_trace_case)
dbstop if error

filenameIndices={'100_0_0_0_0';'0_100_0_0_0';'0_0_100_0_0';'0_0_0_100_0';'0_0_0_0_100'};
load(['SOBO_' defect_case '_FibreType_MinObesrved.mat']);
minObservedValues=[minInitS,minUltStrain,minPdStrain,minUltStrength,minYStrength];
numSensitivityVariables=5;
inputWeights=[0.9 1 1 1 1;1 0.9 1 1 1;1 1 0.9 1 1;1 1 1 0.9 1;1 1 1 1 0.9;1 1 1 1 1;1.1 1 1 1 1;1 1.1 1 1 1;1 1 1.1 1 1;1 1 1 1.1 1;1 1 1 1 1.1];

switch data_trace_case
    case 'VTF_data_trace'
%         load(['MOBO_additive_' defect_case '_results_' filenameIndices{1} '.mat'])
%         XTrace=results.XTrace;
    case {'corrected_VTF_data_trace','corrected_VTF_data_trace_noSmallHybrids'}
        load(['matlab_scripts-second_pass/pristine_VTF_data_corrected.mat'])
%         XTrace=create_table(VTF_data.Carbon_fibre,VTF_data.Glass_fibre,VTF_data.Vc,VTF_data.lf.*1000,VTF_data.SmAvg,VTF_data.G,VTF_data.GiicmAvg);
        XTrace=create_table(Carbon_fibre,Glass_fibre,Vc,lf.*1000,SmAvg,G,GiicmAvg);
%     case {'VTF_data_trace_correctedTest','VTF_data_trace_correctedTest_noSmallHybrids','VTF_data_trace_correctedTest_ignoreSmallHybrids'}
% %         load(['SOBO_single_' defect_case '_results_' filenameIndices{1} '.mat'])
% %         XTrace=results.XTrace;
    otherwise
%         load(['D:/Box Sync/HiPerDuCT PhD/Modelling/Optimisation/results/MOBO/additive/pristine/plots/paretoFronts/' defect_case '_' data_trace_case '_pareto_inputs.mat']) % make sure data_trace_case contains filename you want!
%         load(['matlab_scripts/' defect_case '_observedParetoFront_inputs_alternative.mat']);
%         load(['matlab_scripts/' defect_case '_observedParetoFront_inputs_alternative_corrected.mat']);
%         XTrace=table(Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg); % check lf is consistent
%         MO_case=MO_case;
end
        

for propertyIndex=1:5
    switch data_trace_case
        case {'VTF_data_trace_correctedTest','VTF_data_trace_correctedTest_noSmallHybrids','VTF_data_trace_correctedTest_ignoreSmallHybrids'}
%             load(['SOBO_single_' defect_case '_results_' filenameIndices{propertyIndex} '.mat'])
        otherwise
            load(['results-second_pass/corrected_models/oldMethod_models/corrected/MOBO_additive_' defect_case '_results_' filenameIndices{propertyIndex} '.mat'])
    end
    
%     testin only
%     load(['pristine_VTF_data_corrected_sensitivityInputs.mat'])
%     XTrace=[XTrace;table(Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg)];
    
    if propertyIndex==1
        [meanProperties,stdDevProperties,constraints,stdDevConstraints,constraintViolations]=initialise_data(size(XTrace,1),numSensitivityVariables);
    end   
    
    for sensitivityIndex=1:11
        sensitivityTrace=find_sensitivity_trace(XTrace,inputWeights,sensitivityIndex,data_trace_case);
        if propertyIndex==1
            sensitivityTraceRecord{sensitivityIndex}=sensitivityTrace;
        end
        [meanProperties(:,propertyIndex,sensitivityIndex),stdDevProperties(:,propertyIndex,sensitivityIndex)] = predictObjective(results,sensitivityTrace);
        [constraints(:,propertyIndex,sensitivityIndex),stdDevConstraints(:,propertyIndex,sensitivityIndex)] = predictConstraints(results,sensitivityTrace);
    end
    
    switch data_trace_case % don't need to normalise if using SOBO results
        case {'VTF_data_trace_correctedTest','VTF_data_trace_correctedTest_noSmallHybrids','VTF_data_trace_correctedTest_ignoreSmallHybrids'}
%              meanProperties(:,propertyIndex,:)=meanProperties(:,propertyIndex,:).*-1;
        otherwise
            
            meanProperties(:,propertyIndex,:)=meanProperties(:,propertyIndex,:).*-1.*minObservedValues(propertyIndex);
            stdDevProperties(:,propertyIndex,:)=stdDevProperties(:,propertyIndex,:).*minObservedValues(propertyIndex);
    end
    
end

% StdProperties=meanProperties-1.96.*stdDevProperties; % based on normal distribution with 95% Std
stdProperties=meanProperties-1.*stdDevProperties; 

constraintViolations=reshape(mode(constraints<=0,2).*(-2)+1,size(XTrace,1),2*numSensitivityVariables+1);
constraintViolations=mode(constraintViolations,2); % constraint violations for whole set

meanMeanProperties=mean(meanProperties,3);
stdDevMeanProperties=std(meanProperties,[],3);
CIMeanProperties=meanMeanProperties-2.2281.*stdDevMeanProperties./sqrt(11); % taken from t-table for two-tailed 95% Std with df=n-1=10
stdMeanProperties=meanMeanProperties-1.*stdDevMeanProperties;

[sortMeanProperties,meanSensitivitySortIndex]=sort(meanProperties,3);
minMeanProperties=sortMeanProperties(:,:,1);
minMeanSensitivitySortIndex=meanSensitivitySortIndex(:,:,1);

meanStdProperties=mean(stdProperties,3);
stdDevStdProperties=std(stdProperties,[],3);
% StdStdProperties=meanStdProperties-2.2281.*stdDevStdProperties./sqrt(11); % taken from t-table for two-tailed 95% Std with df=n-1=10
stdStdProperties=meanStdProperties-1.*stdDevStdProperties; % taken from t-table for two-tailed 95% Std with df=n-1=10

[sortStdProperties,stdSensitivitySortIndex]=sort(stdProperties,3);
minStdProperties=sortStdProperties(:,:,1);
minStdSensitivitySortIndex=stdSensitivitySortIndex(:,:,1);

sensitivityCase=findSensitivityCase(meanProperties);

switch data_trace_case
    case 'VTF_data_trace'
%         saveName=['matlab_scripts/local_10pc_' defect_case '_sensitivity_analysis.mat'];
%         save(saveName,'meanProperties','stdDevProperties','stdProperties','constraints','stdDevConstraints','constraintViolations',...
%     'meanMeanProperties','stdDevMeanProperties','stdMeanProperties','meanStdProperties','stdDevStdProperties','stdStdProperties',...
%     'CIMeanProperties','minMeanProperties','minMeanSensitivitySortIndex','minStdProperties','minStdSensitivitySortIndex','sensitivityCase','-v7.3');

    case {'VTF_data_trace_correctedTest','corrected_VTF_data_trace_noSmallHybrids'}
        saveName=['results-second_pass/sensitivityAnalysis/no_small_hybrids/local_10pc_' defect_case '_' data_trace_case '_sensitivity_analysis.mat'];
        save(saveName,'meanProperties','stdDevProperties','stdProperties','constraints','stdDevConstraints','constraintViolations',...
    'meanMeanProperties','stdDevMeanProperties','stdMeanProperties','meanStdProperties','stdDevStdProperties','stdStdProperties',...
    'CIMeanProperties','minMeanProperties','minMeanSensitivitySortIndex','minStdProperties','minStdSensitivitySortIndex','sensitivityCase','-v7.3');

    case {'VTF_data_trace_correctedTest_noSmallHybrids'}
%         saveName=['local_10pc_' defect_case '_sensitivity_analysis_correctedTest_noSmallHybrids.mat'];
%         save(saveName,'meanProperties','stdDevProperties','stdProperties','constraints','stdDevConstraints','constraintViolations',...
%     'meanMeanProperties','stdDevMeanProperties','stdMeanProperties','meanStdProperties','stdDevStdProperties','stdStdProperties',...
%     'CIMeanProperties','minMeanProperties','minMeanSensitivitySortIndex','minStdProperties','minStdSensitivitySortIndex','sensitivityCase','sensitivityTraceRecord','XTrace','-v7.3');


    otherwise
        saveName=['results-second_pass/sensitivityAnalysis/initial_sensitivity_analysis/local_10pc_' defect_case '_' data_trace_case '_sensitivity_analysis.mat'];
        save(saveName,'meanProperties','stdDevProperties','stdProperties','constraints','stdDevConstraints','constraintViolations',...
    'meanMeanProperties','stdDevMeanProperties','stdMeanProperties','meanStdProperties','stdDevStdProperties','stdStdProperties',...
    'CIMeanProperties','minMeanProperties','minMeanSensitivitySortIndex','minStdProperties','minStdSensitivitySortIndex','sensitivityCase','-v7.3');

end
% 
% save(saveName,'meanProperties','stdDevProperties','stdProperties','constraints','stdDevConstraints','constraintViolations',...
%     'meanMeanProperties','stdDevMeanProperties','stdMeanProperties','meanStdProperties','stdDevStdProperties','stdStdProperties',...
%     'CIMeanProperties','minMeanProperties','minMeanSensitivitySortIndex','minStdProperties','minStdSensitivitySortIndex','-v7.3');

end

function sensitivityTrace=find_sensitivity_trace(XTrace,inputWeights,sensitivityIndex,data_trace_case)

Carbon_fibre=categorical(string(table2cell(XTrace(:,'Carbon_fibre'))));
Glass_fibre=categorical(string(table2cell(XTrace(:,'Glass_fibre'))));
% Vc=table2array(XTrace(:,'Vc')).*inputWeights(sensitivityIndex,1);
Vc=table2array(XTrace(:,'Vc'));
switch data_trace_case
    case {'VTF_data_trace_correctedTest','corrected_VTF_data_trace'}
        Vc(Carbon_fibre~=Glass_fibre)=Vc(Carbon_fibre~=Glass_fibre).*inputWeights(sensitivityIndex,1); % only change Vc if Carbon and glass fibres are different
    case {'VTF_data_trace_correctedTest_noSmallHybrids','corrected_VTF_data_trace_noSmallHybrids'}
        Glass_fibre(Vc>=0.95)=Carbon_fibre(Vc>=0.95); % if Vc is greater than 99%, set ADC design to non-hybrid
        Vc(Vc>=0.95)=1;
        Vc(Carbon_fibre~=Glass_fibre)=Vc(Carbon_fibre~=Glass_fibre).*inputWeights(sensitivityIndex,1); % only change Vc if Carbon and glass fibres are different
    case {'VTF_data_trace_correctedTest_ignoreSmallHybrids'}
        Vc(Carbon_fibre~=Glass_fibre & Vc<0.99)=Vc(Carbon_fibre~=Glass_fibre & Vc<0.99).*inputWeights(sensitivityIndex,1); % only change Vc if Carbon and glass fibres are different, and if Vc less than 99%
    otherwise
        Vc=Vc.*inputWeights(sensitivityIndex,1);
end
Vc=min([Vc,ones(size(Vc,1),1)],[],2); % need to make sure Vc never goes over 1!
lf=table2array(XTrace(:,'lf')).*inputWeights(sensitivityIndex,2);
SmAvg=table2array(XTrace(:,'SmAvg')).*inputWeights(sensitivityIndex,3);
G=table2array(XTrace(:,'G')).*inputWeights(sensitivityIndex,4);
GiicmAvg=table2array(XTrace(:,'GiicmAvg')).*inputWeights(sensitivityIndex,5);

sensitivityTrace=table(Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg);

end

function [meanProperties,stdDevProperties,constraints,stdDevConstraints,constraintViolations]=initialise_data(resultsSize,numSensitivityVariables)
    meanProperties=nan(resultsSize,5,numSensitivityVariables*2+1);stdDevProperties=meanProperties;constraints=meanProperties;stdDevConstraints=meanProperties;constraintViolations=nan(resultsSize,numSensitivityVariables*2+1);
end

function sensitivityCase=findSensitivityCase(meanProperties)
sensitivityCase=strings(size(meanProperties,1),size(meanProperties,3));
sensitivityCase(:,1)='low_Vc';
sensitivityCase(:,2)='low_lf';
sensitivityCase(:,3)='low_SmAvg';
sensitivityCase(:,4)='low_G';
sensitivityCase(:,5)='low_GiicmAvg';
sensitivityCase(:,6)='baseline';
sensitivityCase(:,7)='high_Vc';
sensitivityCase(:,8)='high_lf';
sensitivityCase(:,9)='high_SmAvg';
sensitivityCase(:,10)='high_G';
sensitivityCase(:,11)='high_GiicmAvg';
end

function XTrace=create_table(Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg)

 XTrace=table(Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg); % check lf is consistent

end
