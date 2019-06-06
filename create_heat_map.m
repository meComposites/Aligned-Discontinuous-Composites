    function Results = create_heat_map(Results,data,limit_type,zoom_window_height,zoom_window_width,limit_value,plot_title,colormap_limits,expected_value,varargin)
        
%         create more data realisations
%         data=data+flipud(data);
%         data=data+fliplr(data);
%         data=data+rot90(data,1)+rot90(data,2)+rot90(data,3);
%         data=data/16;
        RVE_height = (size(data,1)+1)/2;
        RVE_width = (size(data,2)+1)/2;
        min_zoom_row=RVE_height-zoom_window_height/2;
        max_zoom_row=RVE_height+zoom_window_height/2;
        min_zoom_column=RVE_width-zoom_window_width/2;
        max_zoom_column=RVE_width+zoom_window_width/2;
        
        switch nargin
            case 9
                mean_data=expected_value;
            otherwise
                mean_data=mean(mean(data,'omitnan'),'omitnan');
        end

        
        stdDev_data=sqrt(sum(sum((data-mean_data).^2,'omitnan'),'omitnan')./(sum(sum(1-isnan(data)))-1));

        if stdDev_data == 0
            return
        end
        
        figure()
        switch limit_type
            case 'set limits'
                limits=colormap_limits;
            case 'std dev'
                limits=[mean_data-limit_value*stdDev_data,mean_data+limit_value*stdDev_data];  
            case 'percentage'
                limits=[mean_data-limit_value*mean_data,mean_data+limit_value*mean_data];
            case 'std dev Soraia'
                limits=[mean_data-limit_value,mean_data+limit_value];
        end
        imagesc(data(min_zoom_row:max_zoom_row,min_zoom_column:max_zoom_column),limits)          
        Results.Limits_data=[Results.Limits_data;{{plot_title},{[limits(1,1) mean_data limits(1,2)]}}];
        colormap(redblue)
%         t=title(plot_title);
axis equal tight
axis off
set(gca, 'LooseInset', [0,0,0,0]);
set(gcf,'units','normalized','position',[0.2,0.2,0.25,0.4])
savefig(plot_title)
saveas(gcf,[plot_title '.emf'])       
      
    end