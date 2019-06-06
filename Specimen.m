classdef Specimen
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        RVE % RVE data sets
        Material_properties % material input data
        First_fragment_results % Post-processed results from specimen
        Fracture_fragment_results % Post-processed results from specimen
        Strength_fragment_results % Post-processed results from specimen
        Stress_strain % original stress-strain curve produced by the model (post-processed and averaged across the specimen)
        Stress_strain_fracture % original stress-strain curve but with fracture toughness criterion applied
        Misaligned_stress_strain % Stress_strain_fracture with misalignment applied using equivalent laminate method
        Initial_stiffness % (MPa)
        Ultimate_strain % (%)
        Pseudo_ductile_strain % (%)
        Ultimate_strength %(MPa)
        Yield_strength % 0.1% strain offset yield strength (MPa)
        Critical_RVE_index % RVE index for the RVE who's failure results in catastrophic failure of the specimen
        Critical_cluster % critical cluster resulting in fracture of specimen
        Critical_cluster_loc % rows and columns (1st and 2nd column respectively) for each fibre position in critical cluster
        Critical_cluster_centroid % centroid of the critical cluster (for effects of randomness work)
        Specimen_fracture_strength % the strength of the critical RVE at fracture
        Specimen_fracture_strain % the strain of the critical RVE at fracture
        Fracture_condition % displays 'True' when the specimen failed by fracture, and 'False' otherwise
        Crit_damage_threshold % records the damage threshold at which fracture occurred
        Crit_cluster_size % records number of fibres within critical cluster
        Fragmentation_failures % records the number of fragmentations per fibre
        Debonding_failures % records the number of debonds per fibre
        Softening_failures % records the number of interfaces that soften due to matrix damage
        Fibre_failures % records the number of completely failed fibres 
    end
    
    methods
        
    end
end
