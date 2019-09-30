classdef Fibre_HSR
    %Fibre object with material properties of a fibre
    %   Detailed explanation goes here
    
    properties
        name
        type
        manufacturer
        ef
        EfAvg
        EfCoV
        lfWei
        mfWei
        Sf
        df
    end
    
    methods
        function    obj=Fibre_HSR(name,strain_rate)
            fileID = fopen('Fibre_data.csv');
            fileData = textscan(fileID,'%s %s %s %f %f %f %f %f','delimiter',',');
            fclose(fileID);
            
            if iscategorical(name) ==1
                fibreIndex = find((fileData{2}==name));
            else
                fibreIndex = find(strcmp(fileData{2},name)==1);
            end
            
            obj.name=name;
            obj.type=char(fileData{3}(fibreIndex));
            obj.manufacturer=char(fileData{1}(fibreIndex));
            obj.EfAvg=fileData{4}(fibreIndex)*1000;
            obj.EfCoV=0.05;
            obj.lfWei=fileData{7}(fibreIndex);
            obj.mfWei=fileData{6}(fibreIndex);
            obj.Sf=fileData{5}(fibreIndex);   
            obj.df=fileData{8}(fibreIndex)/1000;
            obj.ef=obj.Sf/(obj.EfAvg)*100;
            obj=determine_Ef_Sf(obj,strain_rate);
        end
        
        function    [obj]=determine_Ef_Sf(obj,strain_rate)
            
            % apply strain rate sensitivity according to {Pawel}
            if (strcmp(obj.type,'E-glass') || strcmp(obj.type,'S-glass')) && isa(strain_rate,'double') && strain_rate>0
                obj.EfAvg=obj.EfAvg+1139.8*strain_rate^0.276; % in MPa
                obj.Sf=obj.Sf+7.71*strain_rate^0.886;
            end
        end
    end
    
end

