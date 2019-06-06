    function create_CDFs(Results,Required_data,Required_axes_titles)
        % loop through each data entry
        for ii = 1:length(Required_data)
            % extract data for each required data entry
            CDF_data=eval(['Results.' Required_data{ii}]);
            % pre-assign cumulative sum to speed-up
            cumul_sum=zeros(1,length(CDF_data));
            % first value s the first value of the data
            cumul_sum(1)=CDF_data(1);
            % cumulative sum
            for jj=2:length(CDF_data)
                cumul_sum(jj)=cumul_sum(jj-1)+CDF_data(jj);
            end
            % normalise to 1
            cumul_sum=(cumul_sum./cumul_sum(end))*(Results.Fracture_count/Results.Specimen_count);
            figure()
            % select axis titles
            switch Required_data{ii}
                case 'Crit_damage_threshold_count'
                    x_axis=0:1/Num_damage_increments:1;
                case 'Crit_cluster_size_count'
                    x_axis=1:RVE_height*RVE_width;
            end
            plot(x_axis,cumul_sum,'-','linewidth',2)
            box off
            set(gca,'ticklabelinterpreter','latex','fontsize',24)
            xlabel(Required_axes_titles{ii},'interpreter','latex','fontsize',24)
            ylabel('CDF','interpreter','latex','fontsize',24)
            set(gcf,'units','normalized','position',[0.2,0.2,0.25,0.4])
            % save figures
            savefig(Required_axes_titles{ii})
            saveas(gcf,[Required_axes_titles{ii} '.emf'])
        end
    end
