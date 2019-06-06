function create_defects_plots(composite_types,defect_types,defect_list,n_repeats)
% Creates defects plots for effects of randomness paper
    property_names={'Ultimate_strain';'Pseudo_ductile_strain';'Ultimate_strength';'Yield_strength';'Initial_stiffness'};
    y_axis_names={'Ultimate strain (\%)';'Pseudo ductile strain (\%)';'Ultimate strength (MPa)';'Yield strength (MPa)';'Initial stiffness (GPa)'};
    defect_names={'Residual fibre fragmentation defects, $\zeta_\mathrm{rff}$ (\%)'; 'Residual interface defects, $\zeta_\mathrm{ri}$ (\%)'; 'Voidage defects, $\zeta_\mathrm{v}$ (\%)'};
    composite_names={'3mm-fibre ADC','3mm-fibre hybrid ADC'};
    colour_palette=[1 164/255 31/255; 0 121/255 212/255; 113/255 174/255 46/255; 138/255 0 0];
    non_present_record=cell(1,1);
    
    data_record=zeros(length(defect_list),n_repeats,length(composite_types),length(property_names),length(defect_types));
    
    for ll = 1:length(composite_types)
        for kk = 1:length(defect_types)
            for jj=1:n_repeats
                for ii=1:length(defect_list)
                    switch defect_types{kk}
                        case 'handling damage'
                            filename=[composite_types{ll} '_' num2str(defect_list(ii)) '_0_0_' num2str(defect_list(ii)) '_0_0_Specimen_' num2str(jj) '.mat'];               
                        case 'interface damage'
                            filename=[composite_types{ll} '_0_' num2str(defect_list(ii)) '_0_0_' num2str(defect_list(ii)) '_0_Specimen_' num2str(jj) '.mat'];
                        case 'voidage'
                            filename=[composite_types{ll} '_0_0_' num2str(defect_list(ii)) '_0_0_' num2str(defect_list(ii)) '_Specimen_' num2str(jj) '.mat'];
                    end
                    try
                        load(filename)
                    catch
                        non_present_record={non_present_record;filename};
                        continue
                    end
                    
                    for mm=1:length(property_names)
                        value=eval(['Specimen_dataset.' property_names{mm}]);
                        if strcmp(property_names{mm},'Yield_strength')==1 && isempty(value)==1
                            value=0;
                        % convert to GPa
                        elseif strcmp(property_names{mm},'Initial_stiffness')==1
                            value=value/1000;
                        end
                        data_record(ii,jj,ll,mm,kk)=value;
                    end
                end
            end
        end
    end
                    
    save('defect_data_record.mat','data_record','-v7.3');
    
    for kk=1:length(defect_types)
        for mm=1:length(property_names)
            figure()
            hold on
            for jj=1:n_repeats
                for ll=1:length(composite_types)
                    plot(defect_list,data_record(:,jj,ll,mm,kk),'LineStyle','none','Marker','x','MarkerSize',8,'LineWidth',2,'MarkerEdgeColor',colour_palette(ll,:))
                end
            end
            xlabel(defect_names{kk},'interpreter','latex','fontsize',18)
            ylabel(y_axis_names{mm},'interpreter','latex','fontsize',18)
            set(gca,'ticklabelinterpreter','latex','fontsize',18)
            h=legend(composite_names);
            set(h,'interpreter','latex','location','best','fontsize',14)
            savefig([defect_types{kk} '_' property_names{mm}]);
        end
    end
end

