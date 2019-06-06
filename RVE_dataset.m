classdef RVE_dataset
    % RVE data set class
    properties
        % variability properties
        RVE_fibre_type % matrix of fibre types (glass or carbon)
        MstocS_array % array of stochastic strengths (all four samples) 
        Fibre_overlap_record % records all initial fibre overlap lengths
 
        
        Matrix_strength_record % matrix strength record around fibre
        
        Inter_fibre_distance_record % inter-fibre distance record
        
        Fibre_stiffness % matrix of fibre stiffness
        Fibre_diameter % matrix of fibre diameters
        
        % results storage
        Stress_strain_unprocessed % original stres-strain curve but without fracture criterion applied 
        Misaligned_stress_strain_unprocessed % misaligned stress strain curve (before fracture)
        Stress_strain_fracture_unprocessed % fracture crierion stress-strain curve stres-strain curve but not processed 
        Fibre_stresses % stresses of each fibre throughout simulation - only use for debugging as it takes up a lot of RAM!
        Fibre_failure_increments % strain increments at which each fibre fails (note, this is not the failure strain, but the model iteration at which it breaks).
        Cumul_crack_count % cumulative crack count for whole RVE
        Carbon_crack_count % cumulative crack count for carbon fibres only
        Glass_crack_count % cumulative crack count for glass fibres only 
        Fibre_crack_count % crack count for each fibre
        Fibre_precrack_count % pre-crack count for each fibre
        Prefailed_interface_count % record of all interface defects that were prefailed at the start of the analysis
        Fibre_partialfailed_count % record of all of the interfaces that have been pre-failed in the RVE (four interface prefailures leads to one fibre prefailure)
        Interface_partialfailed_record % records defective interfaces as above, but specifying which interface has failed
        Fibre_prefailed_count % record of all pre-failed fibres within RVE
        Total_fibre_cracks % total number of fibre cracks from simulation for each fibre
        First_fragment_location % location of first fibre fragmentation location
        Fracture_fragment_location % location of final fibre fragmentation location when fracture criterion is reached - legacy
        Strength_fragment_location % location of final fibre fragmentation location when max strength criterion is reached - legacy
        Xcrit % critical stress at fracture of RVE
        Epscrit % critical strain at fracture of RVE
        Gcrit % critical strain energy release rate calculated at point of fracture
        Clusters_record % a record of all of the clusters within the cross-sectin for each strain increment
        Cluster_size_record % a record of the size of each cluster of broken fibres for every strain increment
        Reserve_factor_record % a record of all the reserve factors for the fracture of every cluster with increasing strain increments
        Cluster_data % size and position of critical cluster within critical RVE [size_rows size_columns row column MinStrainIncrement_Crit] where row and column are the row and column of the top-left corner of the cluster.
        Cumulative_cluster_data % shows how cluster grows and moves around RVE before fracture
        Failed_cluster % shape of failed cluster within RVE
        Nearly_critical_cluster_loc % i and j positions of each fibre within most critical cluster within RVE
        MNumNeighbours_record % matrix of numbers of neighbouring active interfaces for every strain increment
        Mint_record % record of interface stresses for all strain increments
        Fracture_condition % displays 'True' when the specimen failed by fracture, and 'False' otherwise
        Pre_defect_stiffness % the stiffness of the fibres before defects are added
        Damage_state_at_failure % the damage state of all fibres at the critical strain index (at failure)
        Fibre_stagger_distance % the stagger distance of the fibres from the HiPerDiF method
        
        % failure mode identification parameters
        Fragmentation_record % records interface # that led to fragmentation and at what strain
        Failure_mode_record
        
        % parameters from full fracture criterion
        Min_crit_strain_all_fibre_effectiveness % for when fracture toughnes is used when residual stresses are also considered within a cluster. Each min reserve factor saved for every fibre effectiveness, allowing the minimum to be found
        CritStrainClustersAtEffectiveness % shows all of the clusters at the critical strain increment for each fibre effectiveness
        CritClusterAtEffectiveness % shows only critical clusters for each fibre effectiveness increment
        Crit_fibre_effectiveness_index % records the fibre effectiveness index at which fracture occurs
                
        % parameters from FAST fracture method
        Crit_damage_threshold % records the damage threshold at which fracture occurred.
        
        % debugging parameters
        Interface_assert_flag_record
        Interface_assert_data_record
        Fibre_assert_flag_record
        
    end
   
    methods
    end
    
end

