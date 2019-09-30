classdef Fibre
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
        function    obj=Fibre(name)
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
            obj.ef=(obj.Sf/(obj.EfAvg)*100);
        end
    end
    
end

