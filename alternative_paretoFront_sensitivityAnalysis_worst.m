function alternative_paretoFront_sensitivityAnalysis_worst(sensitivityAnalysis_filename,pareto_front_filename,initial_no_small_hybrids_flag,corrected_flag)

load(sensitivityAnalysis_filename);
load(pareto_front_filename);
load(['matlab_scripts-second_pass/pristine_VTF_data_corrected.mat'])
XTrace=create_table(Carbon_fibre,Glass_fibre,Vc,lf.*1000,SmAvg,G,GiicmAvg);

% if sensitivity analysis contains corrected data, get rid of it
switch corrected_flag
    case 'corrected'
%         startingMatfile=matfile('pristine_starting_VTF_data.mat');
%         VcStarting=startingMatfile.Vc;
%         startingSize=size(VcStarting,1);
%         meanProperties=meanProperties(1:size(meanProperties,1)-startingSize,:,:);
end

observedParetoFrontFileIndices=struct('paretoOptimal_ultStrength_initS_fileIndices',paretoOptimal_ultStrength_initS_fileIndices,'paretoOptimal_ultStrength_ultStrain_fileIndices',paretoOptimal_ultStrength_ultStrain_fileIndices,'paretoOptimal_yieldStrength_pdStrain_fileIndices',paretoOptimal_yieldStrength_pdStrain_fileIndices);
fields=fieldnames(observedParetoFrontFileIndices);
[MO_case,fileIndices,XTrace]=find_MO_case(observedParetoFrontFileIndices,fields,XTrace);

meanProperties=meanProperties(fileIndices,:,:);
constraints=constraints(fileIndices,:,:);
sensitivityCase=repmat(["low_Vc","low_lf","low_SmAvg","low_G","low_GiicmAvg","baseline","high_Vc","high_lf","high_SmAvg","high_G","high_GiicmAvg"],size(meanProperties,1),1);

switch initial_no_small_hybrids_flag
    case 'initial'
        save('results-second_pass/sensitivityAnalysis/initial_sensitivity_analysis/paretoFrontSensitivityAnalysis.mat','meanProperties','constraints','MO_case','sensitivityCase','XTrace','-v7.3')
    case 'no_small_hybrids'
        save('results-second_pass/sensitivityAnalysis/no_small_hybrids/paretoFrontSensitivityAnalysis.mat','meanProperties','constraints','MO_case','sensitivityCase','XTrace','-v7.3')
end

end


function XTrace=create_table(Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg)

 XTrace=table(Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg); % check lf is consistent

end


function [MO_case,fileIndices,XTrace]=find_MO_case(observedParetoFrontFileIndices,fields,XTrace)

MO_options=["ultStrength_initS","ultStrength_ultStrain","yieldStrength_pdStrain"];

MO_case=[];
fileIndices=[];
for ii=1:numel(fields)
   MO_case=[MO_case;repmat(MO_options(ii),size(observedParetoFrontFileIndices.(fields{ii}),1),1)];
   fileIndices=[fileIndices;observedParetoFrontFileIndices.(fields{ii})]; 
end

XTrace=XTrace(fileIndices,:);
end

