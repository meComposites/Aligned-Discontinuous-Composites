function gather_VTF_data(defect_case,specific_filename,varargin)
dbstop if error
if nargin ==1
    switch defect_case
        case 'pristine'
%             filedata=dir(['G:/Optimisation/VTF_run_data-second_pass/eConstraintHybrid/' defect_case '/VTF_run_*.mat']);
%             filedata=[filedata;dir(['G:/Optimisation/VTF_run_data-second_pass/' defect_case '/*_run_*.mat'])];
            filedata=dir(['VTF_*run*.mat']);
        case 'defective'
            filedata=dir(['VTF_*run*.mat']);
    end
else
    if strcmp(specific_filename,'variability')==1
        filedata=dir(['G:/Optimisation/VTF_run_data-second_pass/variability/VTF_MO_variability_*.mat']);
    else
        filedata=dir([specific_filename '_run_*.mat']);
    end
end

number_of_files=size(filedata,1);

try
    if nargin == 1
        load([defect_case '_VTF_data.mat'])
    else
        if strcmp(specific_filename,'variability')==1
            load([defect_case '_' specific_filename '_VTF_data.mat'])
        else
            load([specific_filename '_VTF_data.mat'])
        end
    end
    
    Carbon_fibre=string(Carbon_fibre);
    Glass_fibre=string(Glass_fibre);
    VTF_data_size=size(Vc,1);
    data_index=VTF_data_size;
    if VTF_data_size<number_of_files
        Size_increase=number_of_files-VTF_data_size;
        [Filename,UID,Maximise_case,Constraint_case,Slice_number,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,gmAvg,lspecimen,Tu0,Vf,n_rows,n_columns,...
            dc,ec,EcAvg,EcCoV,lcWei,mcWei,Sc,dg,eg,EgAvg,EgCoV,lgWei,mgWei,Sg,Initial_stiffness,Ultimate_strain,...
            Pseudo_ductile_strain,Ultimate_strength,Yield_strength,Max_load_drop,pc_fragmented,pc_vacancy,...
            pg_fragmented,pg_vacancy,p_interface,misalignment_case] = increase_VTF_data_size(...
            Filename,UID,Maximise_case,Constraint_case,Slice_number,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,gmAvg,lspecimen,Tu0,Vf,n_rows,n_columns,...
            dc,ec,EcAvg,EcCoV,lcWei,mcWei,Sc,dg,eg,EgAvg,EgCoV,lgWei,mgWei,Sg,Initial_stiffness,Ultimate_strain,...
            Pseudo_ductile_strain,Ultimate_strength,Yield_strength,Max_load_drop,pc_fragmented,pc_vacancy,...
            pg_fragmented,pg_vacancy,p_interface,misalignment_case,Size_increase);
        
    elseif VTF_data_size>number_of_files
        error('too little files');
    end
    
catch
    disp('no VTF data, starting from scratch!');
    
    [Filename,UID,Maximise_case,Constraint_case,Slice_number,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,gmAvg,lspecimen,Tu0,Vf,n_rows,n_columns,...
        dc,ec,EcAvg,EcCoV,lcWei,mcWei,Sc,dg,eg,EgAvg,EgCoV,lgWei,mgWei,Sg,Initial_stiffness,Ultimate_strain,...
        Pseudo_ductile_strain,Ultimate_strength,Yield_strength,Max_load_drop,pc_fragmented,pc_vacancy,...
        pg_fragmented,pg_vacancy,p_interface,misalignment_case] = initialise_VTF_data(number_of_files);
    
    data_index=0;
    Generation_number=0;
end

for ii=1:number_of_files
    
    isGathered=arrayfun(@(c)strcmp(c,filedata(ii).name),Filename,'UniformOutput',false);
    
    if sum(cell2mat(isGathered))>0
        if sum(cell2mat(isGathered))>1
            error('repeated data');
        end
        continue
    else
        data_index=data_index+1;
        [Filename,UID,Maximise_case,Constraint_case,Slice_number,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,gmAvg,lspecimen,Tu0,Vf,n_rows,n_columns,...
            dc,ec,EcAvg,EcCoV,lcWei,mcWei,Sc,dg,eg,EgAvg,EgCoV,lgWei,mgWei,Sg,Initial_stiffness,Ultimate_strain,...
            Pseudo_ductile_strain,Ultimate_strength,Yield_strength,Max_load_drop,pc_fragmented,pc_vacancy,...
            pg_fragmented,pg_vacancy,p_interface,misalignment_case] = update_VTF_data(...
            Filename,UID,Maximise_case,Constraint_case,Slice_number,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,gmAvg,lspecimen,Tu0,Vf,n_rows,n_columns,...
            dc,ec,EcAvg,EcCoV,lcWei,mcWei,Sc,dg,eg,EgAvg,EgCoV,lgWei,mgWei,Sg,Initial_stiffness,Ultimate_strain,...
            Pseudo_ductile_strain,Ultimate_strength,Yield_strength,Max_load_drop,pc_fragmented,pc_vacancy,...
            pg_fragmented,pg_vacancy,p_interface,misalignment_case,filedata(ii).name,data_index,defect_case);
    end
    
end

Generation_number=[Generation_number;(ones(size(Vc,1)-size(Generation_number,1),1).*(max(Generation_number)+1))];

Carbon_fibre=categorical(Carbon_fibre);
Glass_fibre=categorical(Glass_fibre);

if nargin == 1
    save([defect_case '_VTF_data.mat'],'Filename','UID','Maximise_case','Constraint_case','Slice_number','Generation_number','Carbon_fibre','Glass_fibre','Vc','lf','SmAvg','G','GiicmAvg','gmAvg',...
        'lspecimen','Tu0','Vf','n_rows','n_columns','dc','ec','EcAvg','EcCoV','lcWei','mcWei','Sc',...
        'dg','eg','EgAvg','EgCoV','lgWei','mgWei','Sg','Initial_stiffness','Ultimate_strain',...
        'Pseudo_ductile_strain','Ultimate_strength','Yield_strength','Max_load_drop','pc_fragmented','pc_vacancy',...
        'pg_fragmented','pg_vacancy','p_interface','misalignment_case','-v7.3');
else
    if strcmp(specific_filename,'variability')==1
        save([defect_case '_' specific_filename '_VTF_data.mat'],'Filename','UID','Maximise_case','Constraint_case','Slice_number','Generation_number','Carbon_fibre','Glass_fibre','Vc','lf','SmAvg','G','GiicmAvg','gmAvg',...
            'lspecimen','Tu0','Vf','n_rows','n_columns','dc','ec','EcAvg','EcCoV','lcWei','mcWei','Sc',...
            'dg','eg','EgAvg','EgCoV','lgWei','mgWei','Sg','Initial_stiffness','Ultimate_strain',...
            'Pseudo_ductile_strain','Ultimate_strength','Yield_strength','Max_load_drop','pc_fragmented','pc_vacancy',...
            'pg_fragmented','pg_vacancy','p_interface','misalignment_case','-v7.3');
    else
        save([specific_filename '_VTF_data.mat'],'Filename','UID','Maximise_case','Constraint_case','Slice_number','Generation_number','Carbon_fibre','Glass_fibre','Vc','lf','SmAvg','G','GiicmAvg','gmAvg',...
            'lspecimen','Tu0','Vf','n_rows','n_columns','dc','ec','EcAvg','EcCoV','lcWei','mcWei','Sc',...
            'dg','eg','EgAvg','EgCoV','lgWei','mgWei','Sg','Initial_stiffness','Ultimate_strain',...
            'Pseudo_ductile_strain','Ultimate_strength','Yield_strength','Max_load_drop','pc_fragmented','pc_vacancy',...
            'pg_fragmented','pg_vacancy','p_interface','misalignment_case','-v7.3');
    end
    
end

end

function [Filename,UID,Maximise_case,Constraint_case,Slice_number,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,gmAvg,lspecimen,Tu0,Vf,n_rows,n_columns,...
    dc,ec,EcAvg,EcCoV,lcWei,mcWei,Sc,dg,eg,EgAvg,EgCoV,lgWei,mgWei,Sg,Initial_stiffness,Ultimate_strain,...
    Pseudo_ductile_strain,Ultimate_strength,Yield_strength,Max_load_drop,pc_fragmented,pc_vacancy,...
    pg_fragmented,pg_vacancy,p_interface,misalignment_case] = initialise_VTF_data(number_of_files)
%     Size_increase=3000;

[Filename,UID,Maximise_case,Constraint_case,Carbon_fibre,Glass_fibre]=deal(strings(number_of_files,1));
[Slice_number,Vc,lf,SmAvg,G,GiicmAvg,gmAvg,lspecimen,Tu0,Vf,n_rows,n_columns,dc,ec,EcAvg,EcCoV,lcWei,mcWei,Sc,dg,eg,EgAvg,EgCoV,lgWei,mgWei,Sg,Initial_stiffness,Ultimate_strain,Pseudo_ductile_strain,Ultimate_strength,Yield_strength,Max_load_drop,pc_fragmented,pc_vacancy,pg_fragmented,pg_vacancy,p_interface,misalignment_case]=deal(NaN(number_of_files,1));


% Filename=strings(number_of_files,1);
% UID=strings(number_of_files,1);
% Maximise_case=strings(number_of_files,1);
% Constraint_case=strings(number_of_files,1);
% Slice_number=NaN(number_of_files,1);
% 
% % optimised variables (fibre-type study)
% Carbon_fibre=strings(number_of_files,1);
% Glass_fibre=strings(number_of_files,1);
% Vc=NaN(number_of_files,1);
% lf=NaN(number_of_files,1);
% SmAvg=NaN(number_of_files,1);
% G=NaN(number_of_files,1);
% GiicmAvg=NaN(number_of_files,1);
% 
% % non-optimised / determined properties
% gmAvg=NaN(number_of_files,1);
% lspecimen=NaN(number_of_files,1);
% Tu0=NaN(number_of_files,1);
% Vf=NaN(number_of_files,1);
% n_rows=NaN(number_of_files,1);
% n_columns=NaN(number_of_files,1);
% 
% % individual fibre properties (for blue sky study)
% dc=NaN(number_of_files,1);
% ec=NaN(number_of_files,1);
% EcAvg=NaN(number_of_files,1);
% EcCoV=NaN(number_of_files,1);
% lcWei=NaN(number_of_files,1);
% mcWei=NaN(number_of_files,1);
% Sc=NaN(number_of_files,1);
% 
% dg=NaN(number_of_files,1);
% eg=NaN(number_of_files,1);
% EgAvg=NaN(number_of_files,1);
% EgCoV=NaN(number_of_files,1);
% lgWei=NaN(number_of_files,1);
% mgWei=NaN(number_of_files,1);
% Sg=NaN(number_of_files,1);
% 
% % outputs
% Initial_stiffness=NaN(number_of_files,1);
% Ultimate_strain=NaN(number_of_files,1);
% Pseudo_ductile_strain=NaN(number_of_files,1);
% Ultimate_strength=NaN(number_of_files,1);
% Yield_strength=NaN(number_of_files,1);
% Max_load_drop=NaN(number_of_files,1);
% 
% % defects
% pc_fragmented=NaN(number_of_files,1);
% pc_vacancy=NaN(number_of_files,1);
% pg_fragmented=NaN(number_of_files,1);
% pg_vacancy=NaN(number_of_files,1);
% p_interface=NaN(number_of_files,1);
% misalignment_case=NaN(number_of_files,1);
end

function [Filename,UID,Maximise_case,Constraint_case,Slice_number,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,gmAvg,lspecimen,Tu0,Vf,n_rows,n_columns,...
    dc,ec,EcAvg,EcCoV,lcWei,mcWei,Sc,dg,eg,EgAvg,EgCoV,lgWei,mgWei,Sg,Initial_stiffness,Ultimate_strain,...
    Pseudo_ductile_strain,Ultimate_strength,Yield_strength,Max_load_drop,pc_fragmented,pc_vacancy,...
    pg_fragmented,pg_vacancy,p_interface,misalignment_case] = increase_VTF_data_size(...
    Filename,UID,Maximise_case,Constraint_case,Slice_number,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,gmAvg,lspecimen,Tu0,Vf,n_rows,n_columns,...
    dc,ec,EcAvg,EcCoV,lcWei,mcWei,Sc,dg,eg,EgAvg,EgCoV,lgWei,mgWei,Sg,Initial_stiffness,Ultimate_strain,...
    Pseudo_ductile_strain,Ultimate_strength,Yield_strength,Max_load_drop,pc_fragmented,pc_vacancy,...
    pg_fragmented,pg_vacancy,p_interface,misalignment_case,Size_increase)

%     Size_increase=3000;

Filename=[Filename;strings(Size_increase,1)];
UID=[UID;strings(Size_increase,1)];
Maximise_case=[Maximise_case;strings(Size_increase,1)];
Constraint_case=[Constraint_case;strings(Size_increase,1)];
Slice_number=[Slice_number;NaN(Size_increase,1)];

% optimised variables (fibre-type study)
Carbon_fibre=[Carbon_fibre;strings(Size_increase,1)];
Glass_fibre=[Glass_fibre;strings(Size_increase,1)];
Vc=[Vc;NaN(Size_increase,1)];
lf=[lf;NaN(Size_increase,1)];
SmAvg=[SmAvg;NaN(Size_increase,1)];
G=[G;NaN(Size_increase,1)];
GiicmAvg=[GiicmAvg;NaN(Size_increase,1)];

% non-optimised / determined properties
gmAvg=[gmAvg;NaN(Size_increase,1)];
lspecimen=[lspecimen;NaN(Size_increase,1)];
Tu0=[Tu0;NaN(Size_increase,1)];
Vf=[Vf;NaN(Size_increase,1)];
n_rows=[n_rows;NaN(Size_increase,1)];
n_columns=[n_columns;NaN(Size_increase,1)];

% individual fibre properties (for blue sky study)
dc=[dc;NaN(Size_increase,1)];
ec=[ec;NaN(Size_increase,1)];
EcAvg=[EcAvg;NaN(Size_increase,1)];
EcCoV=[EcCoV;NaN(Size_increase,1)];
lcWei=[lcWei;NaN(Size_increase,1)];
mcWei=[mcWei;NaN(Size_increase,1)];
Sc=[Sc;NaN(Size_increase,1)];

dg=[dg;NaN(Size_increase,1)];
eg=[eg;NaN(Size_increase,1)];
EgAvg=[EgAvg;NaN(Size_increase,1)];
EgCoV=[EgCoV;NaN(Size_increase,1)];
lgWei=[lgWei;NaN(Size_increase,1)];
mgWei=[mgWei;NaN(Size_increase,1)];
Sg=[Sg;NaN(Size_increase,1)];


% outputs
Initial_stiffness=[Initial_stiffness;NaN(Size_increase,1)];
Ultimate_strain=[Ultimate_strain;NaN(Size_increase,1)];
Pseudo_ductile_strain=[Pseudo_ductile_strain;NaN(Size_increase,1)];
Ultimate_strength=[Ultimate_strength;NaN(Size_increase,1)];
Yield_strength=[Yield_strength;NaN(Size_increase,1)];
Max_load_drop=[Max_load_drop;NaN(Size_increase,1)];

% defects
pc_fragmented=[pc_fragmented;NaN(Size_increase,1)];
pc_vacancy=[pc_vacancy;NaN(Size_increase,1)];
pg_fragmented=[pg_fragmented;NaN(Size_increase,1)];
pg_vacancy=[pg_vacancy;NaN(Size_increase,1)];
p_interface=[p_interface;NaN(Size_increase,1)];
misalignment_case=[misalignment_case;NaN(Size_increase,1)];
end

function [Filename,UID,Maximise_case,Constraint_case,Slice_number,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,gmAvg,lspecimen,Tu0,Vf,n_rows,n_columns,...
    dc,ec,EcAvg,EcCoV,lcWei,mcWei,Sc,dg,eg,EgAvg,EgCoV,lgWei,mgWei,Sg,Initial_stiffness,Ultimate_strain,...
    Pseudo_ductile_strain,Ultimate_strength,Yield_strength,Max_load_drop,pc_fragmented,pc_vacancy,...
    pg_fragmented,pg_vacancy,p_interface,misalignment_case] = update_VTF_data(...
    Filename,UID,Maximise_case,Constraint_case,Slice_number,Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg,gmAvg,lspecimen,Tu0,Vf,n_rows,n_columns,...
    dc,ec,EcAvg,EcCoV,lcWei,mcWei,Sc,dg,eg,EgAvg,EgCoV,lgWei,mgWei,Sg,Initial_stiffness,Ultimate_strain,...
    Pseudo_ductile_strain,Ultimate_strength,Yield_strength,Max_load_drop,pc_fragmented,pc_vacancy,...
    pg_fragmented,pg_vacancy,p_interface,misalignment_case,filedata_name,data_index,defect_case)

Filename(data_index)=filedata_name;

UID_start_index=strfind(Filename(data_index),'_run_');
UID_end_index=strfind(Filename(data_index),'.mat');
UID(data_index)=filedata_name(UID_start_index+5:UID_end_index-1);
Defect_case_index=strfind(Filename(data_index),['_' defect_case '_']);
if contains(Filename(data_index),'MOBO')==1
    Maximise_index=strfind(Filename(data_index),'_max_');
    Constraint_index=strfind(Filename(data_index),'_cons_');
    Constraint_slice_string=filedata_name(Constraint_index+6:Defect_case_index-1);
    Slice_index=strfind(Constraint_slice_string,'_'); 
    Maximise_case(data_index)=filedata_name(Maximise_index+5:Constraint_index-1);
    Constraint_case(data_index)=Constraint_slice_string(1:Slice_index(end)-1);
    Slice_number(data_index)=str2double(Constraint_slice_string(Slice_index(end)+1:end));
else
    Maximise_case(data_index)='';
    Constraint_case(data_index)='SOBO';
    Slice_number(data_index)=0;
end


load(filedata_name);

Carbon_fibre(data_index)=Specimen_dataset.Material_properties.Carbon_fibre.name;
Glass_fibre(data_index)=Specimen_dataset.Material_properties.Glass_fibre.name;
Vc(data_index)=Specimen_dataset.Material_properties.Vc;
lf(data_index)=Specimen_dataset.Material_properties.lf;
SmAvg(data_index)=Specimen_dataset.Material_properties.SmAvg;
G(data_index)=Specimen_dataset.Material_properties.G;
GiicmAvg(data_index)=Specimen_dataset.Material_properties.GiicmAvg;

gmAvg(data_index)=Specimen_dataset.Material_properties.gmAvg;
lspecimen(data_index)=Specimen_dataset.Material_properties.lspecimen;
Tu0(data_index)=Specimen_dataset.Material_properties.Tu0;
Vf(data_index)=Specimen_dataset.Material_properties.Vf;
n_rows(data_index)=size(Specimen_dataset.RVE{1}.RVE_fibre_type,1);
n_columns(data_index)=size(Specimen_dataset.RVE{1}.RVE_fibre_type,2);

dc(data_index)=Specimen_dataset.Material_properties.dc;
ec(data_index)=Specimen_dataset.Material_properties.ec;
EcAvg(data_index)=Specimen_dataset.Material_properties.EcAvg;
EcCoV(data_index)=Specimen_dataset.Material_properties.EcCoV;
lcWei(data_index)=Specimen_dataset.Material_properties.lcWei;
mcWei(data_index)=Specimen_dataset.Material_properties.mcWei;
Sc(data_index)=Specimen_dataset.Material_properties.Sc;

dg(data_index)=Specimen_dataset.Material_properties.dg;
eg(data_index)=Specimen_dataset.Material_properties.eg;
EgAvg(data_index)=Specimen_dataset.Material_properties.EgAvg;
EgCoV(data_index)=Specimen_dataset.Material_properties.EgCoV;
lgWei(data_index)=Specimen_dataset.Material_properties.lgWei;
mgWei(data_index)=Specimen_dataset.Material_properties.mgWei;
Sg(data_index)=Specimen_dataset.Material_properties.Sg;

% outputs
Initial_stiffness(data_index)=Specimen_dataset.Initial_stiffness;
Ultimate_strain(data_index)=Specimen_dataset.Ultimate_strain;
Pseudo_ductile_strain(data_index)=Specimen_dataset.Pseudo_ductile_strain;
Ultimate_strength(data_index)=Specimen_dataset.Ultimate_strength;
if isempty(Specimen_dataset.Yield_strength)==1
    Yield_strength(data_index)=NaN;
else
    Yield_strength(data_index)=Specimen_dataset.Yield_strength;
end
if isempty(Specimen_dataset.Max_load_drop)==1
    Max_load_drop(data_index) = determine_max_load_drop(Specimen_dataset.Stress_strain_fracture);
else
    Max_load_drop(data_index)=Specimen_dataset.Max_load_drop;
end

% defects
if isempty(Specimen_dataset.Material_properties.pc_fragmented)==1
    Specimen_dataset.Material_properties.pc_fragmented=0;
    Specimen_dataset.Material_properties.pc_vacancy=0;
    Specimen_dataset.Material_properties.pg_fragmented=0;
    Specimen_dataset.Material_properties.pg_vacancy=0;
    Specimen_dataset.Material_properties.p_interface=0;
    Specimen_dataset.Material_properties.misalignment_case=0;
    save(filedata_name,'Specimen_dataset');
end

pc_fragmented(data_index)=Specimen_dataset.Material_properties.pc_fragmented;
pc_vacancy(data_index)=Specimen_dataset.Material_properties.pc_vacancy;
pg_fragmented(data_index)=Specimen_dataset.Material_properties.pg_fragmented;
pg_vacancy(data_index)=Specimen_dataset.Material_properties.pg_vacancy;
p_interface(data_index)=Specimen_dataset.Material_properties.p_interface;
misalignment_case(data_index)=Specimen_dataset.Material_properties.misalignment_case;
end

function Max_load_drop = determine_max_load_drop(Stress_strain)
% finds the maximum load drop in the final stress-strain curve (part of
% optimisation objective function)
[~,~,~,load_drops]=findpeaks(Stress_strain(:,2));

if isempty(load_drops)==1 || max(load_drops)<=0
    Max_load_drop=0;
else
    Max_load_drop=max(load_drops);
end

end