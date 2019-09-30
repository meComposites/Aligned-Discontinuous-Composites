function [num_files]=gather_MO_pareto_inputs(defect_case,gather_flag,determine_pareto_case)

dbstop if error

% gather data if required
switch gather_flag
    case 'gather'
        % gather data and save into .mat file
        gather_VTF_data(defect_case);
    otherwise
end

switch determine_pareto_case
    case 'determine_paretos'
        % find pareto optimals and associated filenames
        [paretoOptimal_ultStrength_initS,paretoOptimal_ultStrength_ultStrain,paretoOptimal_yieldStrength_pdStrain,paretoOptimal_ultStrength_initS_filenames,paretoOptimal_ultStrength_ultStrain_filenames,paretoOptimal_yieldStrength_pdStrain_filenames]=...
            find_pareto_optimals(defect_case);
    otherwise
        load(['results-second_pass/observedParetoFronts/' determine_pareto_case '.mat']);
end

% concatenate filenames ontop of one another
filenames=[paretoOptimal_ultStrength_initS_filenames;paretoOptimal_ultStrength_ultStrain_filenames;paretoOptimal_yieldStrength_pdStrain_filenames];
num_files=size(filenames,1);
num_ultStrength_initS=size(paretoOptimal_ultStrength_initS_filenames,1);
num_ultStrength_ultStrain=size(paretoOptimal_ultStrength_ultStrain_filenames,1);
num_yieldStrength_pdStrain=size(paretoOptimal_ultStrength_initS_filenames,1);

% initialise variables to speed things up
[Filename,MO_case,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,pc_fragmented,pc_vacancy,pg_fragmented,pg_vacancy,p_interface,misalignment_case,Initial_stiffness,Ultimate_strain,Pseudo_ductile_strain,Ultimate_strength,Yield_strength]=...
    initialise_VTF_data(defect_case,num_files,num_ultStrength_initS,num_ultStrength_ultStrain,num_yieldStrength_pdStrain,paretoOptimal_ultStrength_initS,paretoOptimal_ultStrength_ultStrain,paretoOptimal_yieldStrength_pdStrain);

% get table of input data
data_table=get_table(defect_case);

for ii=1:num_files
    
    isPresent=arrayfun(@(c)strcmp(c,filenames{ii}),data_table.Filename,'UniformOutput',false);
    
    if sum(cell2mat(isPresent))>0
        if sum(cell2mat(isPresent))>1
            error('repeated data');
        end
        
        table_index=find(cell2mat(isPresent));
        
        Filename(ii)=table2array(data_table(table_index,'Filename'));
        Carbon_fibre(ii)=table2array(data_table(table_index,'Carbon_fibre'));
        Glass_fibre(ii)=table2array(data_table(table_index,'Glass_fibre'));
        Vc(ii)=table2array(data_table(table_index,'Vc'));
        lf(ii)=table2array(data_table(table_index,'lf'));
        SmAvg(ii)=table2array(data_table(table_index,'SmAvg'));
        G(ii)=table2array(data_table(table_index,'G'));
        GiicmAvg(ii)=table2array(data_table(table_index,'GiicmAvg'));
        
    else
        
        error(['file ' num2str(ii) 'not found']);
    end
    
    
    
end

Carbon_fibre=categorical(Carbon_fibre);
Glass_fibre=categorical(Glass_fibre);
lf=lf*1000;

switch determine_pareto_case
    case 'determine_paretos'
        save(['matlab_scripts/' defect_case '_MO_pareto_variability_inputs.mat'],'num_files','Filename','MO_case','Carbon_fibre','Glass_fibre','Vc','lf','SmAvg','G','GiicmAvg',...
        'pc_fragmented','pc_vacancy','pg_fragmented','pg_vacancy','p_interface','misalignment_case','Initial_stiffness','Ultimate_strain','Pseudo_ductile_strain','Ultimate_strength','Yield_strength','-v7.3');
    otherwise
        save(['D:/Box Sync/HiPerDuCT PhD\Modelling/Optimisation/results/MOBO/additive/pristine/plots/paretoFronts/' defect_case '_' determine_pareto_case '_pareto_inputs.mat'],'num_files','Filename','MO_case','Carbon_fibre','Glass_fibre','Vc','lf','SmAvg','G','GiicmAvg',...
        'pc_fragmented','pc_vacancy','pg_fragmented','pg_vacancy','p_interface','misalignment_case','Initial_stiffness','Ultimate_strain','Pseudo_ductile_strain','Ultimate_strength','Yield_strength','-v7.3');
end
end

function [Filename,MO_case,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,pc_fragmented,pc_vacancy,pg_fragmented,pg_vacancy,p_interface,misalignment_case,Initial_stiffness,Ultimate_strain,Pseudo_ductile_strain,Ultimate_strength,Yield_strength]=...
    initialise_VTF_data(defect_case,number_of_files,num_ultStrength_initS,num_ultStrength_ultStrain,num_yieldStrength_pdStrain,paretoOptimal_ultStrength_initS,paretoOptimal_ultStrength_ultStrain,paretoOptimal_yieldStrength_pdStrain)

Filename=strings(number_of_files,1);
MO_case=strings(number_of_files,1);

% optimised variables (fibre-type study)
Carbon_fibre=strings(number_of_files,1);
Glass_fibre=strings(number_of_files,1);
Vc=NaN(number_of_files,1);
lf=NaN(number_of_files,1);
SmAvg=NaN(number_of_files,1);
G=NaN(number_of_files,1);
GiicmAvg=NaN(number_of_files,1);

% defects
pc_fragmented=zeros(number_of_files,1).*(strcmp(defect_case,'pristine')==1)+10.*ones(number_of_files,1).*(strcmp(defect_case,'defective')==1);
pc_vacancy=pc_fragmented;
pg_fragmented=pc_fragmented;
pg_vacancy=pc_fragmented;
p_interface=pc_fragmented;
misalignment_case=zeros(number_of_files,1);

% outputs
Initial_stiffness=NaN(size(Vc,1),1);Ultimate_strain=Initial_stiffness;Pseudo_ductile_strain=Initial_stiffness;Ultimate_strength=Initial_stiffness;Yield_strength=Initial_stiffness;

size_initS=size(paretoOptimal_ultStrength_initS,1);
size_ultStrain=size(paretoOptimal_ultStrength_ultStrain,1);

Initial_stiffness(1:size_initS)=paretoOptimal_ultStrength_initS(:,1);
Ultimate_strength(1:size_initS)=paretoOptimal_ultStrength_initS(:,2);
Ultimate_strain(size_initS+1:size_initS+size_ultStrain)=paretoOptimal_ultStrength_ultStrain(:,1);
Ultimate_strength(size_initS+1:size_initS+size_ultStrain)=paretoOptimal_ultStrength_ultStrain(:,2);
Pseudo_ductile_strain(size_initS+size_ultStrain+1:end)=paretoOptimal_yieldStrength_pdStrain(:,1);
Yield_strength(size_initS+size_ultStrain+1:end)=paretoOptimal_yieldStrength_pdStrain(:,2);

% multi-objective optimisation case
MO_case(1:num_ultStrength_initS)="ultStrength_initS";
MO_case(num_ultStrength_initS+1:num_ultStrength_initS+num_ultStrength_ultStrain)="ultStrength_ultStrain";
MO_case(num_ultStrength_initS+num_ultStrength_ultStrain+1:end)="yieldStrength_pdStrain";
end

function data_table=get_table(defect_case)
% load VTF data and store relevant data in table (this avoids the problem
% of variables having the same names)
load([defect_case '_VTF_data.mat'])
data_table=table(Filename,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg);
end
