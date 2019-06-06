classdef Material_properties
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dc % carbon fibre diameter
        ec % ultimate strain of carbon fibre (estimated) (%)
        EcAvg % average carbon fibre stiffness (MPa)
        EcCoV % coefficient of variation for carbon fibre stiffness
        lcWei % carbon fibre Weibull reference length (mm)
        mcWei % carbon fibre Weibull shape parameter
        Sc % carbon fibre Weibull strength scale parameter
        
        dg % glass fibre diameter
        eg % ultimate strain of glass fibre (estimated) (%)
        EgAvg % average glass fibre stiffness (MPa)
        EgCoV % coefficient of variation for glass fibre stiffness 
        lgWei % glass fibre Weibull reference length (mm)
        mgWei % glass fibre Weibull shape parameter
        Sg % glass fibre Weibull strength scale parameter
        
        G % matrix shear stiffness
        GiicmAvg % average matrix mode-II fracture toughness
        gmAvg % average matrix max shear strain (bilinear law)
        SmAvg % average matrix shear strength (bilinear law)
        SmCoV % matrix strength coefficient of variation
        Tu0 % matrix frictional shear stress during fibre pull-out
        
        Vc % carbon ratio (by area or volume)
        Vf % volume fraction (by area or volume)
    end
    
    methods
        function obj = Material_properties(dc,ec,EcAvg,EcCoV,lcWei,mcWei,Sc,dg,eg,EgAvg,EgCoV,lgWei,mgWei,Sg,G,GiicmAvg,gmAvg,SmAvg,SmCoV,Tu0,Vc,Vf)
            
            if nargin < 22
                error('not enough material input data')
            elseif nargin > 22
                error('too much material input data')
            end
            
            obj.dc=dc;
            obj.ec=ec;
            obj.EcAvg=EcAvg;
            obj.EcCoV=EcCoV;
            obj.lcWei=lcWei;
            obj.mcWei=mcWei;
            obj.Sc=Sc;
            
            obj.dg=dg;
            obj.eg=eg;
            obj.EgAvg=EgAvg;
            obj.EgCoV=EgCoV;
            obj.lgWei=lgWei;
            obj.mgWei=mgWei;
            obj.Sg=Sg;
            
            obj.G=G;
            obj.GiicmAvg=GiicmAvg; 
            obj.gmAvg=gmAvg; 
            obj.SmAvg=SmAvg; 
            obj.SmCoV=SmCoV;
            obj.Tu0=Tu0; 
        
            obj.Vc=Vc; 
            obj.Vf=Vf; 
            
        end
    end
    
end

