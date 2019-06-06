classdef Constituent_properties
    % Material properties used for the virtual testing framework
    
    properties
        % matrix
        E11m
        E22m
        G12m
        v12m
        
        % glass fibres
        E11g
        E22g
        G12g
        v12g
        
        % carbon fibres
        E11c
        E22c
        G12c
        v12c
        
        % composite
        Vf
        Vc
        
    end
    
    methods
        function obj = Constituent_properties(E11m,E22m,G12m,v12m,E11g,E22g,G12g,v12g,E11c,E22c,G12c,v12c,Vf,Vc)
                % matrix
        obj.E11m=E11m;
        obj.E22m=E22m;
        obj.G12m=G12m;
        obj.v12m=v12m;
        
        % glass fibres
        obj.E11g=E11g;
        obj.E22g=E22g;
        obj.G12g=G12g;
        obj.v12g=v12g;
        
        % carbon fibres
%         E11c=860000;
        obj.E11c=E11c;
        obj.E22c=E22c;
        obj.G12c=G12c;
        obj.v12c=v12c;
        
        % composite
        obj.Vf=Vf;
        obj.Vc=Vc;
        
        end
        
    end
    
end

