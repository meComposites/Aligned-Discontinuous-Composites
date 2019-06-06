classdef Results_EffectsOfRandomness < handle
    % Class to store post-processed results from experiment
    %   Detailed explanation goes here
    
    properties
        Canvas % the canvas on whih all of the RVEs are superimposed
        Interface_canvas % canvas for nterfaces only (to account for NaNs at edges)
        Specimen_count % the number of specimens that were post-processed
        Fracture_count % the number of specimens that fractured
        
        %% variability properties
        RVE_fibre_type % matrix of fibre types (glass or carbon)
        std_dev_RVE_fibre_type
        MstocS_array % array of stochastic strengths
        MstocS_array_norm % array of stochastic strengths (normalised against expected value for each fibre type)
        MstocS_array_norm_carbon % array of carbon fibre stochastic strengths (normalised against expected value for each fibre type)
        std_dev_MstocS_array_norm_carbon
        MstocS_array_norm_glass % array of glass fibre stochastic strengths (normalised against expected value for each fibre type)
        std_dev_MstocS_array_norm_glass
        MstocS_array_carbon % array of carbon fibre stochastic strengths (not normalised to avoid alphac errors)
        MstocS_array_glass % array of glass fibre stochastic strengths (not normalised to avoid alphag errors)
        
        Overlap_length % all four overlap lengths considered
        std_dev_overlap_length
        
        Matrix_strength % all four interfaces considered
        std_dev_matrix_strength
        
        Inter_fibre_distance % all four interfaces considered
        std_dev_inter_fibre_distance
        
        Fibre_stiffness % matrix of fibre stiffness
        Fibre_stiffness_norm % matrix of fibre stiffness (normalised against expected vlue for each fibre type)
        std_dev_fibre_stiffness_norm
        
        Fibre_diameter % matrix of fibre diameters
        Fibre_diameter_norm % matrix of fibre diameters (normalised against expected value for each fibre type)
        
        %% defects parameters
        Fibre_precrack_count % pre-crack count for each fibre
        std_dev_fibre_precrack_count
        Fibre_partialfailed_count % record of all of the interfaces that have been pre-failed in the RVE (four interface prefailures leads to one fibre prefailure)
        std_dev_fibre_partialfailed_count
        Fibre_prefailed_count % record of all pre-failed fibres within RVE
        std_dev_fibre_prefailed_count
        
        %% CDFs
        Crit_damage_threshold_count % counts the number of times each damage threshold is the critical damage threshold for the specimens
        Crit_cluster_size_count % counts the number of times each cluster is a certain size in the analsis
        
        %% Critical cluster parameters
        Critical_cluster_centroid % heat mp of all of the centroids of the critical clusters
        
        %% Failure modes
        Fragmentation_count % matrix of fibre fragmentations at the fracture strain
        Softening_count % matrix of softening events at the fracture strain
        Debonding_count % matrix of debonds at the fracture strain
        Fibre_failure_count % matrix of fibre failures at the fracture strain
        
        %% plot data
        Std_dev_data % records the calculated standard deviaiton of all of each variability source
        Limits_data % the name of the plot, and the colourmap min, mean, and max values
        
    end
    
    methods
%         function obj = Results_EffectsOfRandomness()
%             obj;
%         end 
%         
%         function Effects_of_randomness(obj,Specimen_array)
%             obj.Fibre_type_distribution = Specimen_array{1,1}.RVE{1,1}.RVE_fibre_type;
%             
%         end
    end
    
end

