classdef Material_properties
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Carbon_fibre % carbon fibre data (for fibre type BO)
        dc % carbon fibre diameter
        ec % ultimate strain of carbon fibre (estimated) (%)
        EcAvg % average carbon fibre stiffness (MPa)
        EcCoV % coefficient of variation for carbon fibre stiffness
        lcWei % carbon fibre Weibull reference length (mm)
        mcWei % carbon fibre Weibull shape parameter
        Sc % carbon fibre Weibull strength scale parameter
        
        Glass_fibre % glass fibre data (for fibre type BO)
        dg % glass fibre diameter
        eg % ultimate strain of glass fibre (estimated) (%)
        EgAvg % average glass fibre stiffness (MPa)
        EgCoV % coefficient of variation for glass fibre stiffness
        lgWei % glass fibre Weibull reference length (mm)
        mgWei % glass fibre Weibull shape parameter
        Sg % glass fibre Weibull strength scale parameter
        
        G % matrix shear stiffnes
        GiicmAvg % average matrix mode-II fracture toughness
        gmAvg % average matrix max shear strain (bilinear law)
        SmAvg % average matrix shear strength (bilinear law)
        SmCoV % matrix strength coefficient of variation
        Tu0 % matrix frictional shear stress during fibre pull-out
        
        Vc % carbon ratio (by area or volume)
        Vf % volume fraction (by area or volume)
        lf % fibre lengths
        lspecimen % specimen length
        
        pc_fragmented % fraction of pre-fragmented carbon fibres
        pc_vacancy % fraction of carbon fibres used for vacancy defects
        pg_fragmented % fraction of pre-fragmented glass fibres
        pg_vacancy % fraction of glass fibres used for vacancy defects
        p_interface % fraction of pre-failed interfaces
        misalignment_case % type of misalignment data applied (0= no misalignment, 1=HiPerDiF{Yu2014}, 2={Sanadi})

        Gcrit_c % carbon fibre translaminar fracture toughness
        Gcrit_g % glass fibre translaminar fracture toughness
    end
    
    methods
        function obj = Material_properties(Carbon_fibre,dc,ec,EcAvg,EcCoV,...
                lcWei,mcWei,Sc,Glass_fibre,dg,eg,EgAvg,EgCoV,lgWei,mgWei,...
                Sg,G,GiicmAvg,gmAvg,SmAvg,SmCoV,Tu0,Vc,Vf,lf,lspecimen,...
                pc_fragmented,pc_vacancy,pg_fragmented,pg_vacancy,p_interface,...
                misalignment_case,Gcrit_c,Gcrit_g)
            
            obj.Carbon_fibre=Carbon_fibre;
            obj.dc=dc;
            obj.ec=ec;
            obj.EcAvg=EcAvg;
            obj.EcCoV=EcCoV;
            obj.lcWei=lcWei;
            obj.mcWei=mcWei;
            obj.Sc=Sc;
            
            obj.Glass_fibre=Glass_fibre;
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
            obj.lf=lf;
            obj.lspecimen=lspecimen;
            
            obj.pc_fragmented=pc_fragmented;
            obj.pc_vacancy=pc_vacancy;
            obj.pg_fragmented=pg_fragmented;
            obj.pg_vacancy=pg_vacancy;
            obj.p_interface=p_interface; 
            obj.misalignment_case=misalignment_case;
            
            obj.Gcrit_c=Gcrit_c;
            obj.Gcrit_g=Gcrit_g;
            
        end
    end
    
end

