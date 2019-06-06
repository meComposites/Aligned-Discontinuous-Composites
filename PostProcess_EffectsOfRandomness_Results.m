function [Results] = PostProcess_EffectsOfRandomness_Results(n_Specimens,filename,Composite_type,specimen_list,varargin)
% FJ - a function to amalgamate all of the results from the HPC to give
% contour plots and statistical data about the specimens tested

% if ~contains(Composite_type,'long')==1
%     lf = 0.5;
% else
%     lf = 3;
% end

% load([Composite_type '_inputs.mat']);

Num_damage_increments=100;

switch Composite_type
    case 'HSC'
        Material=Material_properties(0.007,1.9307,225000,0.05,25,5.7119,4344,0.01,1.9307,225000,0.05,25,5.7119,4344,1300,0.8,0.04,60,0.15,10,1,0.55);
    case 'HMC_EG'
        Material=Material_properties(0.007,0.3988,860000,0.05,1,5,3430,0.007,3.2877,73000,0.05,25,5.7119,2400,1500,0.8,0.04,60,0.15,10,49/149,0.3); 
    case 'HSC_Joel'
        Material=Material_properties(0.01,1.9307,225000,0.05,25,5.7119,4344,0.01,1.9307,225000,0.05,25,5.7119,4344,1000,0.8,0.04,100,0.15,10,1,0.36);
end

dbstop if error

lf = input('Please input lf: ');
lspecimen=lf;

X_avg_c=Material.Sc*gamma(1+1/Material.mcWei)*((Material.lcWei/(lf/4))^(1/Material.mcWei));
X_avg_g=Material.Sg*gamma(1+1/Material.mgWei)*((Material.lgWei/(lf/4))^(1/Material.mgWei));

% load the first specimen dataset to get the specimen size (note: all
% specimens must be the same size for this to work correctly)
data=[];
for start_index = 1:n_Specimens+1
    try
        data=load([filename '_Specimen_' num2str(start_index) '.mat']);
    catch
        if start_index > n_Specimens
            error('cannot find data')
        end
    end
    if isempty(data)~=1
        break
    end
end
    

% create global variables for manipulation of data
RVE_height=size(data.Specimen_dataset.RVE{1}.RVE_fibre_type,1);
RVE_width=size(data.Specimen_dataset.RVE{1}.RVE_fibre_type,2);
Canvas_height=2*RVE_height-1;
Canvas_width=2*RVE_width-1;
% initialise the results arrays
[Results] = Initialise_results_arrays();

% loop through all specimen data
for n=1:n_Specimens
    % load specimen data
    try 
        if nargin == 4 && find((specimen_list==n),1,'first')==0
            continue
        else 
            data=load([filename '_Specimen_' num2str(n) '.mat']);
            
        end
    catch 
        continue
    end
    
    % update results as you go
    [Results] = Update_results_arrays(Results,data.Specimen_dataset);
    Results.Specimen_count=Results.Specimen_count+1;
    % clear data to save RAM
    clear data centroid
end

% divide results by number of times each specimen hd a fibre in that
% section of the canvas
[Results] = Finalise_results_arrays(Results);

% plot results
Plot_results(Results);

% save results in .mat file
save([filename '_FinalResults_EffectsOfRandomness.mat'],'Results','-v7.3')


% Initialise results arays with NaNs - NaNs are used to prevent glass fibres
% being mistaken for zeros etc when the results matrices are
% amalgamated
    function [Results] = Initialise_results_arrays()
        % declare the results strct
        Results = Results_EffectsOfRandomness;
        Results.Specimen_count=0;
        Results.Fracture_count=0;
        Initial_array=nan(Canvas_height,Canvas_width);
        
        Results.Canvas=Initial_array;
        Results.Interface_canvas=Initial_array;
        
        %variability
                
        Results.std_dev_RVE_fibre_type=0;        
        Results.std_dev_MstocS_array_norm_carbon=0;
        Results.std_dev_MstocS_array_norm_glass=0;
        Results.std_dev_overlap_length=0;
        Results.std_dev_matrix_strength=0;
        Results.std_dev_inter_fibre_distance=0;
        Results.std_dev_fibre_stiffness_norm=0;
        Results.std_dev_fibre_precrack_count=0;
        Results.std_dev_fibre_partialfailed_count=0;
        Results.std_dev_fibre_prefailed_count=0;
        
        Results.RVE_fibre_type=Initial_array;
%         Results.MstocS_array=Initial_array;
%         Results.MstocS_array_norm=Initial_array;
        Results.MstocS_array_norm_carbon=Initial_array;
        Results.MstocS_array_norm_glass=Initial_array;
%         Results.MstocS_array_carbon=Initial_array;
%         Results.MstocS_array_glass=Initial_array;
        
        Results.Overlap_length=Initial_array;
        
        Results.Matrix_strength=Initial_array;
        
        Results.Inter_fibre_distance=Initial_array;
        
%         Results.Fibre_stiffness=Initial_array;
        Results.Fibre_stiffness_norm=Initial_array;
%         Results.Fibre_diameter=Initial_array;
        
        % defects
        Results.Fibre_precrack_count=Initial_array; 
        Results.Fibre_partialfailed_count=Initial_array; 
        Results.Fibre_prefailed_count=Initial_array;
        
        % CDFs
        Results.Crit_damage_threshold_count=zeros(1,Num_damage_increments+1);
        Results.Crit_cluster_size_count=zeros(1,RVE_height*RVE_width);
%         Results.Crit_damage_threshold_count=-1*ones(1,n_Specimens);
%         Results.Crit_cluster_size_count=-1*ones(1,n_Specimens);

        % cluster centroid heat map
        Results.Critical_cluster_centroid=zeros(RVE_height,RVE_width);
        
        % failure modes
        Results.Fragmentation_count=Initial_array;
        Results.Softening_count=Initial_array;
        Results.Debonding_count=Initial_array;
        Results.Fibre_failure_count=Initial_array;
    end


    function [Results] = Update_results_arrays(Results,Specimen_dataset)
        
        % extract data from datasets to make things easier
        Critical_RVE_index=Specimen_dataset.Critical_RVE_index;
        RVE_data=Specimen_dataset.RVE{Critical_RVE_index};
        
                % Cluster size
        Crit_cluster_size=sum(sum(Specimen_dataset.Critical_cluster));
        if Crit_cluster_size==0
%             Crit_cluster_size=RVE_height*RVE_width;
%             warning('cluster size may not be defined properly');
            return
        else 
            if isempty(Specimen_dataset.Critical_cluster_centroid)==1
                Specimen_dataset.Critical_cluster_centroid(1,1)=round(sum(Specimen_dataset.Critical_cluster_loc(:,1),1)/size(Specimen_dataset.Critical_cluster_loc(:,1),1)+rand()-0.5);
                Specimen_dataset.Critical_cluster_centroid(1,2)=round(sum(Specimen_dataset.Critical_cluster_loc(:,2),1)/size(Specimen_dataset.Critical_cluster_loc(:,2),1)+rand()-0.5);
            end
            Results.Fracture_count=Results.Fracture_count+1;
            i_centroid=Specimen_dataset.Critical_cluster_centroid(1);
            j_centroid=Specimen_dataset.Critical_cluster_centroid(2);
        end
        
        % pre-calculate co-ordinates of results updates to same time and to
        % keep things clear
        i_min=RVE_height-i_centroid+1;
        i_max=2*RVE_height-i_centroid;
        j_min=RVE_width-j_centroid+1;
        j_max=2*RVE_width-j_centroid;
        
        if i_centroid<1 || j_centroid<1 || i_centroid>RVE_height || j_centroid>RVE_width
            warning('critical cluster lies outside of RVE boundary')
            return
        end
        
        NaN_temp=ones(RVE_height,RVE_width)*NaN;
        Fibre_temp=RVE_data.RVE_fibre_type;
        relative_strength_carbon=NaN_temp;
        relative_strength_glass=NaN_temp;
        MstocS_array=RVE_data.MstocS_array(:,:,1);
        
        % Use cumulative sum to add the new dataset results to the old ones
        % (but only take last slice,as this is the most up-to-date). cumsum
        % used because you can omit NaNs in summations, making the mean
        % values correct.
        temp=cumsum(cat(3,Results.Canvas(i_min:i_max,j_min:j_max),ones(RVE_height,RVE_width)),3,'omitnan');
        Results.Canvas(i_min:i_max,j_min:j_max)=temp(:,:,end);
        % same for interface canvas, but account for NaNs at edges
        % top interface
        temp=cumsum(cat(3,Results.Interface_canvas(i_min:i_max,j_min:j_max),[NaN(1,RVE_width);ones(RVE_height-1,RVE_width)]),3,'omitnan');
        Results.Interface_canvas(i_min:i_max,j_min:j_max)=temp(:,:,end);
        % right interface
        temp=cumsum(cat(3,Results.Interface_canvas(i_min:i_max,j_min:j_max),[ones(RVE_height,RVE_width-1),NaN(RVE_height,1)]),3,'omitnan');
        Results.Interface_canvas(i_min:i_max,j_min:j_max)=temp(:,:,end);
        % bottom interface
        temp=cumsum(cat(3,Results.Interface_canvas(i_min:i_max,j_min:j_max),[ones(RVE_height-1,RVE_width);NaN(1,RVE_width)]),3,'omitnan');
        Results.Interface_canvas(i_min:i_max,j_min:j_max)=temp(:,:,end);
        % left interface
        temp=cumsum(cat(3,Results.Interface_canvas(i_min:i_max,j_min:j_max),[NaN(RVE_height,1),ones(RVE_height,RVE_width-1)]),3,'omitnan');
        Results.Interface_canvas(i_min:i_max,j_min:j_max)=temp(:,:,end);
        
        % variability
        mean_RVE_fibre_type=mean(mean(RVE_data.RVE_fibre_type,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.RVE_fibre_type(i_min:i_max,j_min:j_max),RVE_data.RVE_fibre_type),3,'omitnan');
        Results.RVE_fibre_type(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_RVE_fibre_type = Results.std_dev_RVE_fibre_type+sum(sum((RVE_data.RVE_fibre_type-mean_RVE_fibre_type).^2,'omitnan'),'omitnan');
        
%         temp=cumsum(cat(3,Results.MstocS_array(i_min:i_max,j_min:j_max),RVE_data.MstocS_array(:,:,1)),3,'omitnan');
%         Results.MstocS_array(i_min:i_max,j_min:j_max)=temp(:,:,end);
        % MstocS normalised
        relative_strength = bsxfun(@rdivide,RVE_data.MstocS_array(:,:,1),((RVE_data.RVE_fibre_type==1).*X_avg_c)+(((1-RVE_data.RVE_fibre_type)==1).*X_avg_g));
%         temp=cumsum(cat(3,Results.MstocS_array_norm(i_min:i_max,j_min:j_max),relative_strength),3,'omitnan');
%         Results.MstocS_array_norm(i_min:i_max,j_min:j_max)=temp(:,:,end);
        % carbon fibres only
%         Canvas_index=(find(Fibre_temp)+((j_min-1)*Canvas_height+i_min-1));
        
        Canvas_index=Inf(Canvas_height,Canvas_width);
        Canvas_index(i_min:i_max,j_min:j_max)=Fibre_temp;
        Canvas_index=find(Canvas_index==1);
        mean_MstocS_array_norm_carbon=mean(mean(relative_strength(Fibre_temp==1),'omitnan'),'omitnan');
        temp=cumsum(cat(2,Results.MstocS_array_norm_carbon(Canvas_index),relative_strength(Fibre_temp==1)),2,'omitnan');
        Results.MstocS_array_norm_carbon(Canvas_index)=temp(:,2);
        Results.std_dev_MstocS_array_norm_carbon=Results.std_dev_MstocS_array_norm_carbon+sum(sum((relative_strength(Fibre_temp==1)/mean_MstocS_array_norm_carbon-1).^2,'omitnan'),'omitnan');
        
        % glass fibres only
        Canvas_index=Inf(Canvas_height,Canvas_width);
        Canvas_index(i_min:i_max,j_min:j_max)=Fibre_temp;
        Canvas_index=find(Canvas_index==0);
        mean_MstocS_array_norm_glass=mean(mean(relative_strength(Fibre_temp==0),'omitnan'),'omitnan');
        temp=cumsum(cat(2,Results.MstocS_array_norm_glass(Canvas_index),relative_strength(Fibre_temp==0)),2,'omitnan');
        Results.MstocS_array_norm_glass(Canvas_index)=temp(:,2);
        Results.std_dev_MstocS_array_norm_glass=Results.std_dev_MstocS_array_norm_glass+sum(sum((relative_strength(Fibre_temp==0)/mean_MstocS_array_norm_glass-1).^2,'omitnan'),'omitnan');
        
%         % non-normalised fibre strengths (to avoid alphac and alphag
%         % problems)
%         % carbon fibres only
%         Canvas_index=Inf(Canvas_height,Canvas_width);
%         Canvas_index(i_min:i_max,j_min:j_max)=Fibre_temp;
%         Canvas_index=find(Canvas_index==1);
%         temp=cumsum(cat(2,Results.MstocS_array_norm_carbon(Canvas_index),MstocS_array(Fibre_temp==1)),2,'omitnan');
%         Results.MstocS_array_carbon(Canvas_index)=temp(:,2);
%         
%         % glass fibres only
%         Canvas_index=Inf(Canvas_height,Canvas_width);
%         Canvas_index(i_min:i_max,j_min:j_max)=Fibre_temp;
%         Canvas_index=find(Canvas_index==0);
%         temp=cumsum(cat(2,Results.MstocS_array_norm_glass(Canvas_index),MstocS_array(Fibre_temp==0)),2,'omitnan');
%         Results.MstocS_array_glass(Canvas_index)=temp(:,2);
        
        % add all four contributions when looking at interface properties
        temp=cumsum(cat(3,Results.Overlap_length(i_min:i_max,j_min:j_max),RVE_data.Fibre_overlap_record(:,:,1)),3,'omitnan');
        Results.Overlap_length(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_overlap_length=Results.std_dev_overlap_length+sum(sum((RVE_data.Fibre_overlap_record(:,:,1)-lf/4).^2,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Overlap_length(i_min:i_max,j_min:j_max),RVE_data.Fibre_overlap_record(:,:,2)),3,'omitnan');
        Results.Overlap_length(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_overlap_length=Results.std_dev_overlap_length+sum(sum((RVE_data.Fibre_overlap_record(:,:,2)-lf/4).^2,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Overlap_length(i_min:i_max,j_min:j_max),RVE_data.Fibre_overlap_record(:,:,3)),3,'omitnan');
        Results.Overlap_length(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_overlap_length=Results.std_dev_overlap_length+sum(sum((RVE_data.Fibre_overlap_record(:,:,3)-lf/4).^2,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Overlap_length(i_min:i_max,j_min:j_max),RVE_data.Fibre_overlap_record(:,:,4)),3,'omitnan');
        Results.Overlap_length(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_overlap_length=Results.std_dev_overlap_length+sum(sum((RVE_data.Fibre_overlap_record(:,:,4)-lf/4).^2,'omitnan'),'omitnan');
        
        temp=cumsum(cat(3,Results.Matrix_strength(i_min:i_max,j_min:j_max),RVE_data.Matrix_strength_record(:,:,1)),3,'omitnan');
        Results.Matrix_strength(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_matrix_strength=Results.std_dev_matrix_strength+sum(sum((RVE_data.Matrix_strength_record(:,:,1)-Material.SmAvg).^2,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Matrix_strength(i_min:i_max,j_min:j_max),RVE_data.Matrix_strength_record(:,:,2)),3,'omitnan');
        Results.Matrix_strength(i_min:i_max,j_min:j_max)=temp(:,:,end);        
        Results.std_dev_matrix_strength=Results.std_dev_matrix_strength+sum(sum((RVE_data.Matrix_strength_record(:,:,2)-Material.SmAvg).^2,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Matrix_strength(i_min:i_max,j_min:j_max),RVE_data.Matrix_strength_record(:,:,3)),3,'omitnan');
        Results.Matrix_strength(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_matrix_strength=Results.std_dev_matrix_strength+sum(sum((RVE_data.Matrix_strength_record(:,:,3)-Material.SmAvg).^2,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Matrix_strength(i_min:i_max,j_min:j_max),RVE_data.Matrix_strength_record(:,:,4)),3,'omitnan');
        Results.Matrix_strength(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_matrix_strength=Results.std_dev_matrix_strength+sum(sum((RVE_data.Matrix_strength_record(:,:,4)-Material.SmAvg).^2,'omitnan'),'omitnan');
        
        mean_inter_fibre_distance=mean(mean(mean(RVE_data.Inter_fibre_distance_record,'omitnan'),'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Inter_fibre_distance(i_min:i_max,j_min:j_max),RVE_data.Inter_fibre_distance_record(:,:,1)),3,'omitnan');
        Results.Inter_fibre_distance(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_inter_fibre_distance=Results.std_dev_inter_fibre_distance+sum(sum((RVE_data.Inter_fibre_distance_record(:,:,1)-mean_inter_fibre_distance).^2,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Inter_fibre_distance(i_min:i_max,j_min:j_max),RVE_data.Inter_fibre_distance_record(:,:,2)),3,'omitnan');
        Results.Inter_fibre_distance(i_min:i_max,j_min:j_max)=temp(:,:,end);        
        Results.std_dev_inter_fibre_distance=Results.std_dev_inter_fibre_distance+sum(sum((RVE_data.Inter_fibre_distance_record(:,:,2)-mean_inter_fibre_distance).^2,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Inter_fibre_distance(i_min:i_max,j_min:j_max),RVE_data.Inter_fibre_distance_record(:,:,3)),3,'omitnan');
        Results.Inter_fibre_distance(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_inter_fibre_distance=Results.std_dev_inter_fibre_distance+sum(sum((RVE_data.Inter_fibre_distance_record(:,:,3)-mean_inter_fibre_distance).^2,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Inter_fibre_distance(i_min:i_max,j_min:j_max),RVE_data.Inter_fibre_distance_record(:,:,4)),3,'omitnan');
        Results.Inter_fibre_distance(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_inter_fibre_distance=Results.std_dev_inter_fibre_distance+sum(sum((RVE_data.Inter_fibre_distance_record(:,:,4)-mean_inter_fibre_distance).^2,'omitnan'),'omitnan');
        
%         temp=cumsum(cat(3,Results.Fibre_stiffness(i_min:i_max,j_min:j_max),RVE_data.Fibre_stiffness),3,'omitnan');
%         Results.Fibre_stiffness(i_min:i_max,j_min:j_max)=temp(:,:,end);
        % stiffness normalised
        relative_stiffness = bsxfun(@rdivide,RVE_data.Fibre_stiffness(:,:,1),((RVE_data.RVE_fibre_type==1).*Material.EcAvg)+(((1-RVE_data.RVE_fibre_type)==1).*Material.EgAvg));
        temp=cumsum(cat(3,Results.Fibre_stiffness_norm(i_min:i_max,j_min:j_max),relative_stiffness),3,'omitnan');
        Results.Fibre_stiffness_norm(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_fibre_stiffness_norm=Results.std_dev_fibre_stiffness_norm+sum(sum((relative_stiffness-1).^2,'omitnan'),'omitnan');
        
%         temp=cumsum(cat(3,Results.Fibre_diameter(i_min:i_max,j_min:j_max),RVE_data.Fibre_diameter),3,'omitnan');
%         Results.Fibre_diameter(i_min:i_max,j_min:j_max)=temp(:,:,end);
        
        % defects
        mean_fibre_precrack_count=mean(mean(RVE_data.Fibre_precrack_count,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Fibre_precrack_count(i_min:i_max,j_min:j_max),RVE_data.Fibre_precrack_count),3,'omitnan');
        Results.Fibre_precrack_count(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_fibre_precrack_count=Results.std_dev_fibre_precrack_count+sum(sum((RVE_data.Fibre_precrack_count-mean_fibre_precrack_count).^2,'omitnan'),'omitnan');
        
        mean_fibre_partialfailed_count=mean(mean(mean(RVE_data.Interface_partialfailed_record,'omitnan'),'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Fibre_partialfailed_count(i_min:i_max,j_min:j_max),RVE_data.Interface_partialfailed_record(:,:,1)),3,'omitnan');
        Results.Fibre_partialfailed_count(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_fibre_partialfailed_count=Results.std_dev_fibre_partialfailed_count+sum(sum((RVE_data.Interface_partialfailed_record(:,:,1)-mean_fibre_partialfailed_count).^2,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Fibre_partialfailed_count(i_min:i_max,j_min:j_max),RVE_data.Interface_partialfailed_record(:,:,2)),3,'omitnan');
        Results.Fibre_partialfailed_count(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_fibre_partialfailed_count=Results.std_dev_fibre_partialfailed_count+sum(sum((RVE_data.Interface_partialfailed_record(:,:,2)-mean_fibre_partialfailed_count).^2,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Fibre_partialfailed_count(i_min:i_max,j_min:j_max),RVE_data.Interface_partialfailed_record(:,:,3)),3,'omitnan');
        Results.Fibre_partialfailed_count(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_fibre_partialfailed_count=Results.std_dev_fibre_partialfailed_count+sum(sum((RVE_data.Interface_partialfailed_record(:,:,3)-mean_fibre_partialfailed_count).^2,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Fibre_partialfailed_count(i_min:i_max,j_min:j_max),RVE_data.Interface_partialfailed_record(:,:,4)),3,'omitnan');
        Results.Fibre_partialfailed_count(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_fibre_partialfailed_count=Results.std_dev_fibre_partialfailed_count+sum(sum((RVE_data.Interface_partialfailed_record(:,:,4)-mean_fibre_partialfailed_count).^2,'omitnan'),'omitnan');
        
        mean_fibre_prefailed_count=mean(mean(RVE_data.Fibre_prefailed_count,'omitnan'),'omitnan');
        temp=cumsum(cat(3,Results.Fibre_prefailed_count(i_min:i_max,j_min:j_max),RVE_data.Fibre_prefailed_count),3,'omitnan');
        Results.Fibre_prefailed_count(i_min:i_max,j_min:j_max)=temp(:,:,end);
        Results.std_dev_fibre_prefailed_count=Results.std_dev_fibre_prefailed_count+sum(sum((RVE_data.Fibre_prefailed_count-mean_fibre_prefailed_count).^2,'omitnan'),'omitnan');

        % CDFs
        % Critical damage threshold
        Crit_damage_threshold=Specimen_dataset.RVE{Specimen_dataset.Critical_RVE_index}.Crit_damage_threshold;
        Results.Crit_damage_threshold_count(round(Crit_damage_threshold*Num_damage_increments+1))=Results.Crit_damage_threshold_count(round(Crit_damage_threshold*Num_damage_increments+1))+1;
        
        Results.Crit_cluster_size_count(Crit_cluster_size)=Results.Crit_cluster_size_count(Crit_cluster_size)+1;
        
        %% heat maps
        Results.Critical_cluster_centroid=Results.Critical_cluster_centroid+(Specimen_dataset.Critical_cluster./Crit_cluster_size);
        
        %% failure modes
        temp=cumsum(cat(3,Results.Fragmentation_count(i_min:i_max,j_min:j_max),Specimen_dataset.Fragmentation_failures),3,'omitnan');
        Results.Fragmentation_count(i_min:i_max,j_min:j_max)=temp(:,:,end);
        temp=cumsum(cat(3,Results.Softening_count(i_min:i_max,j_min:j_max),Specimen_dataset.Softening_failures),3,'omitnan');
        Results.Softening_count(i_min:i_max,j_min:j_max)=temp(:,:,end);
        temp=cumsum(cat(3,Results.Debonding_count(i_min:i_max,j_min:j_max),Specimen_dataset.Debonding_failures),3,'omitnan');
        Results.Debonding_count(i_min:i_max,j_min:j_max)=temp(:,:,end);
        temp=cumsum(cat(3,Results.Fibre_failure_count(i_min:i_max,j_min:j_max),Specimen_dataset.Fibre_failures),3,'omitnan');
        Results.Fibre_failure_count(i_min:i_max,j_min:j_max)=temp(:,:,end);
    end

    function [Results] = Finalise_results_arrays(Results)       
        
        Results.Canvas=Results.Canvas;
%         Results.Std_dev_data=[{'Property name'},{'Std. dev.'}];
        
        % variability
        Results.RVE_fibre_type=bsxfun(@rdivide,Results.RVE_fibre_type,Results.Canvas);
        Results.std_dev_RVE_fibre_type=sqrt(Results.std_dev_RVE_fibre_type/sum(sum(Results.Canvas,'omitnan'),'omitnan'));
%         Results.MstocS_array=bsxfun(@rdivide,Results.MstocS_array,Results.Canvas);
%         Results.MstocS_array_norm=bsxfun(@rdivide,Results.MstocS_array_norm,Results.Canvas);
%         Results.MstocS_array_norm_carbon=bsxfun(@rdivide,Results.MstocS_array_norm_carbon,Results.Canvas);
%         Results.MstocS_array_norm_glass=bsxfun(@rdivide,Results.MstocS_array_norm_glass,Results.Canvas);
        Results.MstocS_array_norm_carbon=bsxfun(@rdivide,Results.MstocS_array_norm_carbon,bsxfun(@times,Results.Canvas,Results.RVE_fibre_type));
        Results.MstocS_array_norm_carbon=bsxfun(@rdivide,Results.MstocS_array_norm_carbon,mean(mean(Results.MstocS_array_norm_carbon,'omitnan'),'omitnan'));
        Results.std_dev_MstocS_array_norm_carbon=sqrt(Results.std_dev_MstocS_array_norm_carbon/sum(sum(bsxfun(@times,Results.Canvas,Results.RVE_fibre_type),'omitnan'),'omitnan'));
        
        Results.MstocS_array_norm_glass=bsxfun(@rdivide,Results.MstocS_array_norm_glass,bsxfun(@times,Results.Canvas,(1-Results.RVE_fibre_type)));
        Results.MstocS_array_norm_glass=bsxfun(@rdivide,Results.MstocS_array_norm_glass,mean(mean(Results.MstocS_array_norm_glass,'omitnan'),'omitnan'));
        Results.std_dev_MstocS_array_norm_glass=sqrt(Results.std_dev_MstocS_array_norm_glass/sum(sum(bsxfun(@times,Results.Canvas,(1-Results.RVE_fibre_type)),'omitnan'),'omitnan'));
%         Results.MstocS_array_carbon=bsxfun(@rdivide,Results.MstocS_array_norm_carbon,bsxfun(@times,Results.Canvas,Results.RVE_fibre_type));
%         Results.MstocS_array_glass=bsxfun(@rdivide,Results.MstocS_array_norm_glass,bsxfun(@times,Results.Canvas,(1-Results.RVE_fibre_type)));
        
        % divide by 4 times the number of canvases for interface
        % variability properties because four interfaces
        Results.Overlap_length=bsxfun(@rdivide,Results.Overlap_length,Results.Interface_canvas);
        Results.std_dev_overlap_length=sqrt(Results.std_dev_overlap_length/sum(sum(Results.Interface_canvas,'omitnan'),'omitnan'));

        Results.Matrix_strength=bsxfun(@rdivide,Results.Matrix_strength,Results.Interface_canvas);
        Results.std_dev_matrix_strength=sqrt(Results.std_dev_matrix_strength/sum(sum(Results.Interface_canvas,'omitnan'),'omitnan'));
        
        Results.Inter_fibre_distance=bsxfun(@rdivide,Results.Inter_fibre_distance,Results.Interface_canvas);
        Results.std_dev_inter_fibre_distance=sqrt(Results.std_dev_inter_fibre_distance/sum(sum(Results.Interface_canvas,'omitnan'),'omitnan'));
        
%         Results.Fibre_stiffness=bsxfun(@rdivide,Results.Fibre_stiffness,Results.Canvas);
        Results.Fibre_stiffness_norm=bsxfun(@rdivide,Results.Fibre_stiffness_norm,Results.Canvas);
        Results.std_dev_fibre_stiffness_norm=sqrt(Results.std_dev_fibre_stiffness_norm/sum(sum(Results.Canvas,'omitnan'),'omitnan'));
%         Results.Fibre_diameter=bsxfun(@rdivide,Results.Fibre_diameter,Results.Canvas);
        
        % defects
        Results.Fibre_precrack_count=bsxfun(@rdivide,Results.Fibre_precrack_count,Results.Canvas);
        Results.std_dev_fibre_precrack_count=sqrt(Results.std_dev_fibre_precrack_count/sum(sum(Results.Canvas,'omitnan'),'omitnan'));
        
        Results.Fibre_partialfailed_count=bsxfun(@rdivide,Results.Fibre_partialfailed_count,Results.Interface_canvas);
        Results.std_dev_fibre_partialfailed_count=sqrt(Results.std_dev_fibre_partialfailed_count/sum(sum(Results.Interface_canvas,'omitnan'),'omitnan'));
        
        Results.Fibre_prefailed_count=bsxfun(@rdivide,Results.Fibre_prefailed_count,Results.Canvas);
        Results.std_dev_fibre_prefailed_count=sqrt(Results.std_dev_fibre_prefailed_count/sum(sum(Results.Canvas,'omitnan'),'omitnan'));
        
        Results.Critical_cluster_centroid=Results.Critical_cluster_centroid./max(max(Results.Canvas));
        
        % failure modes
        Results.Fragmentation_count=bsxfun(@rdivide,Results.Fragmentation_count,Results.Canvas);
        % divide by two because interface failures are double counted
        Results.Softening_count=bsxfun(@rdivide,Results.Softening_count,2*Results.Canvas);
        Results.Debonding_count=bsxfun(@rdivide,Results.Debonding_count,2*Results.Canvas);
        Results.Fibre_failure_count=bsxfun(@rdivide,Results.Fibre_failure_count,Results.Canvas);

    end

    function Plot_results(Results)
        Results.Limits_data=[];
        Results.Limits_data=cell(1,2);
        Results.Limits_data{1,1}='Plot title';
        Results.Limits_data{1,2}='Colourmap limits';
        
        figure()
        t=title('Canvas coverage');
        imagesc(Results.Canvas)
        axis equal tight
        
        % plot reslts - basic atm but can be improved
        A=Results.Canvas;
        A=A+flipud(A);
        A=A+fliplr(A);
        A=A+rot90(A,1)+rot90(A,2)+rot90(A,3);
        figure()
        t=title('Canvas manpulation');
        imagesc(A)
        axis equal tight
        %         axis off
        
        zoom_window_height = 20;
        zoom_window_width = 20;
        
        % heat maps - only show up if the std dev of the data isn't zero
        % (i.e. defects plots ony show if there are defects in the
        % specimen)
        % variability
%         Results=create_heat_map(Results,Results.RVE_fibre_type,'std dev',zoom_window_height,zoom_window_width,3,'RVE fibre type');
%         Results=create_heat_map(Results,Results.MstocS_array_norm,'std dev',zoom_window_height,zoom_window_width,3,'Fibre stochastic strength (normalised)');
%         Results=create_heat_map(Results,Results.MstocS_array_norm_carbon,'std dev',zoom_window_height,zoom_window_width,3,'Carbon fibre stochastic strength (normalised)');
%         Results=create_heat_map(Results,Results.MstocS_array_norm_glass,'std dev',zoom_window_height,zoom_window_width,3,'Glass fibre stochastic strength (normalised)');
%         Results=create_heat_map(Results,Results.Overlap_length,'std dev',zoom_window_height,zoom_window_width,3,'Overlap length');
%         Results=create_heat_map(Results,Results.Matrix_strength,'std dev',zoom_window_height,zoom_window_width,3,'Matrix strength');
%         Results=create_heat_map(Results,Results.Inter_fibre_distance,'std dev',zoom_window_height,zoom_window_width,3,'Inter-fibre distance');
%         Results=create_heat_map(Results,Results.Fibre_stiffness_norm,'std dev',zoom_window_height,zoom_window_width,3,'Fibre stiffness (normalised)');
        
        Results=create_heat_map(Results,Results.RVE_fibre_type,'std dev Soraia',zoom_window_height,zoom_window_width,Results.std_dev_RVE_fibre_type,'RVE fibre type');
        Results=create_heat_map(Results,Results.MstocS_array_norm_carbon,'std dev Soraia',zoom_window_height,zoom_window_width,Results.std_dev_MstocS_array_norm_carbon,'Carbon fibre stochastic strength (normalised)');
        Results=create_heat_map(Results,Results.MstocS_array_norm_glass,'std dev Soraia',zoom_window_height,zoom_window_width,Results.std_dev_MstocS_array_norm_glass,'Glass fibre stochastic strength (normalised)');
        Results=create_heat_map(Results,Results.Overlap_length,'std dev Soraia',zoom_window_height,zoom_window_width,Results.std_dev_overlap_length,'Overlap length');
        Results=create_heat_map(Results,Results.Matrix_strength,'std dev Soraia',zoom_window_height,zoom_window_width,Results.std_dev_matrix_strength,'Matrix strength');
        Results=create_heat_map(Results,Results.Inter_fibre_distance,'std dev Soraia',zoom_window_height,zoom_window_width,Results.std_dev_inter_fibre_distance,'Inter-fibre distance');
        Results=create_heat_map(Results,Results.Fibre_stiffness_norm,'std dev Soraia',zoom_window_height,zoom_window_width,Results.std_dev_fibre_stiffness_norm,'Fibre stiffness (normalised)');

        % defects
        Results=create_heat_map(Results,Results.Fibre_precrack_count,'std dev Soraia',zoom_window_height,zoom_window_width,Results.std_dev_fibre_precrack_count,'Handling damage defects');
        Results=create_heat_map(Results,Results.Fibre_partialfailed_count,'std dev Soraia',zoom_window_height,zoom_window_width,Results.std_dev_fibre_partialfailed_count,'Interface damage defects');
        Results=create_heat_map(Results,Results.Fibre_prefailed_count,'std dev Soraia',zoom_window_height,zoom_window_width,Results.std_dev_fibre_prefailed_count,'Voidage damage defects');
     
        % CDFs for looking at trends in critical damage threshold etc.
%         Results.Crit_damage_threshold_count=Results.Crit_damage_threshold_count(Results.Crit_damage_threshold_count>=0);
%         cdfplot(Results.Crit_damage_threshold_count);
%         Results.Crit_cluster_size_count=Results.Crit_cluster_size_count(Results.Crit_cluster_size_count>=0);
%         cdfplot(Results.Crit_cluster_size_count);
        create_CDFs(Results,{'Crit_damage_threshold_count','Crit_cluster_size_count'},{'Critical damage threshold','Critical cluster size (number of fibres)'})
        
        % heat map for critical cluster centroid
        figure()        
        imagesc(Results.Critical_cluster_centroid,[-1E-03,1E-03])
        colormap(redblue)
        axis equal tight
        axis off
        set(gca, 'LooseInset', [0,0,0,0]);
        set(gcf,'units','normalized','position',[0.2,0.2,0.25,0.4])
        savefig('Critical cluster centroid')
        saveas(gcf,['Critical cluster centroid' '.emf'])
        
        Results=create_heat_map(Results,Results.Fragmentation_count,'set limits',zoom_window_height,zoom_window_width,3,'Fragmentation count',[-4,4]);
        Results=create_heat_map(Results,Results.Softening_count,'set limits',zoom_window_height,zoom_window_width,3,'Softening count',[-2,2]);
        Results=create_heat_map(Results,Results.Debonding_count,'set limits',zoom_window_height,zoom_window_width,3,'Debonding count',[-2,2]);
        Results=create_heat_map(Results,Results.Fibre_failure_count,'set limits',zoom_window_height,zoom_window_width,3,'Fibre failure count',[-1,1]);
    end

    % create CDF pots for each required data entry
    function create_CDFs(Results,Required_data,Required_axes_titles)
        % loop through each data entry
        for ii = 1:length(Required_data)
            % extract data for each required data entry
            CDF_data=eval(['Results.' Required_data{ii}]);
            % pre-assign cumulative sum to speed-up
            cumul_sum=zeros(1,length(CDF_data));
            % first value s the first value of the data
            cumul_sum(1)=CDF_data(1);
            % cumulative sum
            for jj=2:length(CDF_data)
                cumul_sum(jj)=cumul_sum(jj-1)+CDF_data(jj);
            end
            % normalise to 1
            cumul_sum=(cumul_sum./cumul_sum(end))*(Results.Fracture_count/Results.Specimen_count);
            figure()
            % select axis titles
            switch Required_data{ii}
                case 'Crit_damage_threshold_count'
                    x_axis=0:1/Num_damage_increments:1;
                case 'Crit_cluster_size_count'
                    x_axis=1:RVE_height*RVE_width;
            end
            plot(x_axis,cumul_sum,'-','linewidth',2)
            box off
            set(gca,'ticklabelinterpreter','latex','fontsize',24)
            xlabel(Required_axes_titles{ii},'interpreter','latex','fontsize',24)
            ylabel('CDF','interpreter','latex','fontsize',24)
            set(gcf,'units','normalized','position',[0.2,0.2,0.25,0.4])
            % save figures
            savefig(Required_axes_titles{ii})
            saveas(gcf,[Required_axes_titles{ii} '.emf'])
        end
    end


end