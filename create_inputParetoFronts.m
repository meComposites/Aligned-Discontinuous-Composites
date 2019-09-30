%% create Pareto front plots that use the inputs to colour the poits of the Pareto fronts

function create_inputParetoFronts(defect_case,Pareto_type,save_flag,non_hybrid_label_flag)
dbstop if error
switch defect_case
    case 'pristine'
        
        switch Pareto_type
            case 'estimatedMean'
                load('matlab_scripts-second_pass/VTF_data_corrected_struct.mat')
                [Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg]=get_struct_inputs(VTF_data);
                clearvars -except Carbon_fibre Glass_fibre Vc lf SmAvg G GiicmAvg defect_case Pareto_type save_flag non_hybrid_label_flag
                load('results-second_pass/estimatedParetoFronts/estimatedParetoFronts2.mat')
            
            case 'noSmallHybrids'
                load('D:\Box Sync\HiPerDuCT PhD\Modelling\Optimisation\results\MOBO\additive\pristine\plots\sensitvityAnalysis-correctedTest\4 - run sensitivity analysis\local_10pc_pristine_sensitivity_analysis_correctedTest_noSmallHybrids.mat')
                XTrace=modify_XTrace(XTrace,Pareto_type);
                Carbon_fibre=categorical(table2array(XTrace(:,'Carbon_fibre')));
                Glass_fibre=categorical(table2array(XTrace(:,'Glass_fibre')));
                clearvars -except Carbon_fibre Glass_fibre XTrace defect_case Pareto_type save_flag
                load('D:\Box Sync\HiPerDuCT PhD\Modelling\Optimisation\results\MOBO\additive\pristine\plots\sensitvityAnalysis-correctedTest\2 - determine observed Pareto front\estimatedParetoFronts_correctedTest_noSmallHybrids.mat')
        end
end

Material_system = find_material_system(Carbon_fibre,Glass_fibre,non_hybrid_label_flag);
paretoOptimal_ultStrength_initS(:,1)=paretoOptimal_ultStrength_initS(:,1)./1000; % scale to GPa
pareto_front_data=struct('paretoOptimal_ultStrength_initS',paretoOptimal_ultStrength_initS,'paretoOptimal_ultStrength_ultStrain',paretoOptimal_ultStrength_ultStrain,'paretoOptimal_yieldStrength_pdStrain',paretoOptimal_yieldStrength_pdStrain);
data_fields=fieldnames(pareto_front_data);
pareto_front_fileIndices=struct('paretoOptimal_ultStrength_initS_fileIndices',paretoOptimal_ultStrength_initS_fileIndices,'paretoOptimal_ultStrength_ultStrain_fileIndices',paretoOptimal_ultStrength_ultStrain_fileIndices,'paretoOptimal_yieldStrength_pdStrain_fileIndices',paretoOptimal_yieldStrength_pdStrain_fileIndices);
indices_fields=fieldnames(pareto_front_fileIndices);
point_colour_data=struct('Material_system',Material_system,'Vc',Vc,'lf',lf,'SmAvg',SmAvg,'G',G,'GiicmAvg',GiicmAvg);
colour_names=fieldnames(point_colour_data);
colourmap_limits=[0,0;0 1;0.5 12;40 100;1000 1800;0.6 1.0];
colourbar_labels={'Material system','$V_\mathrm{c}$','$l_\mathrm{f}$ $\mathrm{\left(mm\right)}$','$\textbf{E} \left(S_\mathrm{m}\right)$ $\mathrm{\left(MPa\right)}$','$G_\mathrm{m}$ $\mathrm{\left(MPa\right)}$','$\mathcal{G}_\mathrm{II_m}$ $\mathrm{\left(kJ/m^{2}\right)}$'};

x_axis_index=[1,2,3];
y_axis_index=[4,4,5];

xlimits=[0 600;0 5;0 5];

for pareto_front_index=1:length(data_fields)
    % put data n correct poistions in array (makes sure the axes are
    % consistent when using colour_scatter function)
    plot_data=NaN(size(pareto_front_data.(data_fields{pareto_front_index}),1),5);
    plot_data(:,x_axis_index(pareto_front_index))=pareto_front_data.(data_fields{pareto_front_index})(:,1);
    plot_data(:,y_axis_index(pareto_front_index))=pareto_front_data.(data_fields{pareto_front_index})(:,2);
    
    for input_data_index=1:length(colour_names)
        if input_data_index == 1   
            input_data=cellstr(Material_system(pareto_front_fileIndices.(indices_fields{pareto_front_index}))); % use string for array referencing, but convert to cell array for compatibility with colour plot function?
            create_colour_scatter_plot('none',defect_case,'colour_scatter',plot_data,input_data,'matSys',[],x_axis_index(pareto_front_index),y_axis_index(pareto_front_index),'colorcube')
            fig=gcf;
        else
            input_data=point_colour_data.(colour_names{input_data_index})(pareto_front_fileIndices.(indices_fields{pareto_front_index})); % selects the appropriate data (Vc, lf, etc.) for the colour of the plots
            create_colour_scatter_plot('none',defect_case,'colour_scatter',plot_data,input_data,point_colour_data.(colour_names{input_data_index}),[],x_axis_index(pareto_front_index),y_axis_index(pareto_front_index),'redblue',colourmap_limits(input_data_index,:),colourbar_labels{input_data_index})
        end
        
        ylim([0 2000])
        xlim(xlimits(pareto_front_index,:))
        
        switch save_flag
            case 'save'
                plot_filepath=['results-second_pass/inputParetoFronts/' data_fields{pareto_front_index} '/'];
                plot_name=colour_names{input_data_index};
                try
                    savefig([plot_filepath,plot_name]);
                    saveas(gcf,[plot_filepath plot_name '.emf'])
                catch
                    mkdir(plot_filepath);
                    savefig([plot_filepath,plot_name]);
                    saveas(gcf,[plot_filepath plot_name '.emf'])
                end
        end
            
        figure(fig)
    end
end



end

function [Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg] = get_inputs(XTrace)

Carbon_fibre=categorical(table2array(XTrace(:,'Carbon_fibre')));
Glass_fibre=categorical(table2array(XTrace(:,'Glass_fibre')));
Vc=table2array(XTrace(:,'Vc'));
lf=table2array(XTrace(:,'lf'))./1000;
SmAvg=table2array(XTrace(:,'SmAvg'));
G=table2array(XTrace(:,'G'));
GiicmAvg=table2array(XTrace(:,'GiicmAvg'));
                
end

function [Carbon_fibre,Glass_fibre,Vc,lf,SmAvg,G,GiicmAvg] = get_struct_inputs(VTF_data)

Carbon_fibre=VTF_data.Carbon_fibre;
Glass_fibre=VTF_data.Glass_fibre;
Vc=VTF_data.Vc;
lf=VTF_data.lf;
SmAvg=VTF_data.SmAvg;
G=VTF_data.G;
GiicmAvg=VTF_data.GiicmAvg;
                
end

function XTrace=modify_XTrace(XTrace,Pareto_front)

switch Pareto_front
    case 'noSmallHybrids'
        % ASSUMPTION THAT DATA IS ALREADY CORRECTED
        
        % modify XTrace to apply no small hybrids condition (anythng with
        % Vc>0.99 set to non-hybrid
        XTrace(table2array(XTrace(:,'Vc'))>0.99,'Glass_fibre')=XTrace(table2array(XTrace(:,'Vc'))>0.99,'Carbon_fibre');
        Vc=table2array(XTrace(:,'Vc'));
        Vc(Vc>0.99)=1;
        XTrace(:,'Vc')=table(Vc);
end
end

function Material_system = find_material_system(Carbon_fibre,Glass_fibre,non_hybrid_label_flag)

Carbon_fibre= string(Carbon_fibre);
Glass_fibre= string(Glass_fibre);

[Carbon_fibre_type,Glass_fibre_type,Material_system]=deal(cell(length(Carbon_fibre),1));

for ii=1:length(Carbon_fibre)
    if isempty(Carbon_fibre_type{ii})==1
        if max(arrayfun(@(n) strcmp(n, char(Carbon_fibre{ii})), ["K13D","XN-05","XN-90","P120J","P75S"]))==1
            Carbon_fibre_type{ii}='Pitch-CF';
        elseif max(arrayfun(@(n) strcmp(n, char(Carbon_fibre{ii})), ["C100","GF"]))==1
            Carbon_fibre_type{ii}='E-Glass';
        elseif max(arrayfun(@(n) strcmp(n, char(Carbon_fibre{ii})), ["FliteStrand_S_ZT"]))==1
            Carbon_fibre_type{ii}='S-Glass';
        else
            Carbon_fibre_type{ii}='PAN-CF';
        end
    end
    
    if isempty(Glass_fibre_type{ii})==1
        if max(arrayfun(@(n) strcmp(n, char(Glass_fibre{ii})), ["K13D","XN-05","XN-90","P120J","P75S"]))==1
            Glass_fibre_type{ii}='Pitch-CF';
        elseif max(arrayfun(@(n) strcmp(n, char(Glass_fibre{ii})), ["C100","GF"]))==1
            Glass_fibre_type{ii}='E-Glass';
        elseif max(arrayfun(@(n) strcmp(n, char(Glass_fibre{ii})), ["FliteStrand_S_ZT"]))==1
            Glass_fibre_type{ii}='S-Glass';
        else
            Glass_fibre_type{ii}='PAN-CF';
        end
    end
    switch non_hybrid_label_flag
        case 'non-hybrid'
            if strcmp(Carbon_fibre{ii},Glass_fibre{ii})==1
                Material_system{ii} = strjoin(["non-hybrid " Carbon_fibre_type{ii}],"");
            else
                Material_system{ii} = strjoin([Carbon_fibre_type{ii} ", " Glass_fibre_type{ii}],"");
            end
        otherwise    
            Material_system{ii} = strjoin([Carbon_fibre_type{ii} ", " Glass_fibre_type{ii}],"");
    end
end

Material_system=string(Material_system);

end


