function create_colour_scatter_plot(save_flag,defect_type,plot_type,plot_data,input_data,input_name,k_clusters,x_axis_index,y_axis_index,colormap_name,colourmap_limits,colourbar_label,kMeans_data,varargin)

colour_map=eval(['colormap(' colormap_name ')']);
plot_names={'initS','ultStrain','pdStrain','ultStrength','yieldStrength'};
axis_labels={'Initial stiffness (GPa)','Ultimate strain(\%)','Pseudo-ductile strain(\%)','Ultimate strength (MPa)','Yield strength (MPa)'};

if strcmp(input_name,'matSys')~=1
    switch plot_type
        case '3d_scatter'
            clusters=ones(size(plot_data,1),1);
            k_clusters=1;
            colorvalues=[0 121/255 212/255];
        case '3d_kMeans_scatter'
            clusters=kmeans(kMeans_data,k_clusters,'replicates',15);
            interval=round(length(colour_map)/k_clusters);
            colorvalues=colour_map(1:interval:end,:);
        case 'kMeans'
            clusters=kmeans(kMeans_data,k_clusters,'replicates',15);
            interval=round(length(colour_map)/k_clusters);
            colorvalues=colour_map(1:interval:end,:);
        case 'no_clusters'
            clusters=ones(size(plot_data,1),1);
            k_clusters=1;
            colorvalues=[0 121/255 212/255];
        otherwise
            if size(input_data,2)~=1
                error('input_data cannot have more than 1 dimension for this plot type')
            end
            
            k_clusters=size(colour_map,1);
            colorvalues=colour_map;
            
            if nargin>=11
                if isempty(colourmap_limits)==1
                    colourmap_limits=[min(input_data),max(input_data)];
                elseif size(colourmap_limits)~=[1,2]
                    error('colourmap limits should be a 1x2 vector');
                end
                input_interval=(colourmap_limits(2)-colourmap_limits(1))/k_clusters;
                clusters=floor((input_data-colourmap_limits(1))/input_interval);
            else
                input_interval=(max(input_data)-min(input_data))/k_clusters;
                clusters=floor((input_data-min(input_data))/input_interval);
            end
            
            clusters(clusters==0)=1;
            
    end
else
    colour_map=colormap(colorcube);
    %     colour_map=colour_map(randperm(size(colour_map,1)),:);
    input_data=categorical(input_data);
    %     matSysCategories=categories(input_data);
    matSysCategories=create_matSysSategories(); % changed this to have all categories (consistent legend?)
    k_clusters=length(matSysCategories);
    clusters=zeros(length(input_data),1);
    
    for kk = 1:k_clusters
        clusters(:,kk)=arrayfun(@(n) strcmp(char(n), char(matSysCategories(kk)))*kk, input_data);
    end
    
    clusters=max(clusters,[],2);
    
    interval=round(length(colour_map)/k_clusters);
    colorvalues=colour_map(1:interval:end,:);
    
    colour_order=[7 14 11 6 5 15 9 1 13 17 16 18 2 10 19 8 3 12 4 20];
    colorvalues = colorvalues(colour_order(unique(clusters)),:);
    %     matSysCategories=cellstr(string(matSysCategories(colour_order(unique(clusters)))));
    k_clusters=size(unique(clusters),1);
    cluster_indices=unique(clusters);
    matSysCategories=cellstr(string(matSysCategories(cluster_indices)));
    matSysCategories=flipud(matSysCategories); % ONLY WHEN INCREMENTING NUMBER OFCLUSTERS IN DESCENDING ORDER
end


for jj=1:length(x_axis_index)
    fig=figure();
    hold on
    
    if contains(plot_type,'3d')==1
        for ii=1:k_clusters
            %             scatter3(plot_data(clusters==ii,x_axis_index(jj)),plot_data(clusters==ii,y_axis_index(jj)),input_data(clusters==ii),'o','markerfacecolor',colorvalues(ii,:),'markeredgecolor',[0,0,0],'markersize',4);
            scatter3(plot_data(clusters==ii,x_axis_index(jj)),input_data(clusters==ii),plot_data(clusters==ii,y_axis_index(jj)),10,'o','markerfacecolor',colorvalues(ii,:),'markeredgecolor',[0,0,0]);
            view(3)
            xlabel(axis_labels{x_axis_index(jj)},'FontSize',14,'interpreter','latex');
            zlabel(axis_labels{y_axis_index(jj)},'FontSize',14,'interpreter','latex');
            ylabel(colourbar_label,'FontSize',14,'interpreter','latex');
            
        end
    else
        if strcmp(input_name,'matSys')==1
            for ii=length(cluster_indices):-1:1
                %         for ii=k_clusters:-1:1
                %             plot(plot_data(clusters==ii,x_axis_index(jj)),plot_data(clusters==ii,y_axis_index(jj)),'o','markerfacecolor',colorvalues(ii,:),'markeredgecolor',[0,0,0],'markersize',4)
                %                 plot3(plot_data(clusters==cluster_indices(ii),x_axis_index(jj)),plot_data(clusters==cluster_indices(ii),y_axis_index(jj)),rand(size(plot_data(clusters==cluster_indices(ii),y_axis_index(jj)))),'o','markerfacecolor',colorvalues((ii),:),'markeredgecolor',[0,0,0],'markersize',6)
                %                 plot3(plot_data(clusters==cluster_indices(ii),x_axis_index(jj)),plot_data(clusters==cluster_indices(ii),y_axis_index(jj)),zeros(size(plot_data(clusters==cluster_indices(ii),y_axis_index(jj)))),'o','markerfacecolor',colorvalues((ii),:),'markeredgecolor','none','markersize',6)
                plot3(plot_data(clusters==cluster_indices(ii),x_axis_index(jj)),plot_data(clusters==cluster_indices(ii),y_axis_index(jj)),rand(size(plot_data(clusters==cluster_indices(ii),y_axis_index(jj)))),'o','markerfacecolor',colorvalues((ii),:),'markeredgecolor',[0,0,0],'markersize',6)
                %         plot(plot_data(clusters==k_clusters-ii+1,x_axis_index(jj)),plot_data(clusters==k_clusters-ii+1,y_axis_index(jj)),'o','markerfacecolor',colorvalues(k_clusters-ii+1,:),'markeredgecolor',[0,0,0],'markersize',4)
            end
            %             plot(plot_data(:,x_axis_index(jj)),plot_data(:,y_axis_index(jj)),'--','lineWidth',0.5,'color',[127,127,127]./255);
            plot3(plot_data(:,x_axis_index(jj)),plot_data(:,y_axis_index(jj)),ones(size(plot_data(:,y_axis_index(jj))))*5,'--','lineWidth',0.5,'color',[127,127,127]./255);
        else
            for ii=1:k_clusters
                %         for ii=k_clusters:-1:1
                %             plot(plot_data(clusters==ii,x_axis_index(jj)),plot_data(clusters==ii,y_axis_index(jj)),'o','markerfacecolor',colorvalues(ii,:),'markeredgecolor',[0,0,0],'markersize',4)
                plot3(plot_data(clusters==ii,x_axis_index(jj)),plot_data(clusters==ii,y_axis_index(jj)),rand(size(plot_data(clusters==ii,y_axis_index(jj)))),'o','markerfacecolor',colorvalues(ii,:),'markeredgecolor',[0,0,0],'markersize',6)
                %         plot(plot_data(clusters==k_clusters-ii+1,x_axis_index(jj)),plot_data(clusters==k_clusters-ii+1,y_axis_index(jj)),'o','markerfacecolor',colorvalues(k_clusters-ii+1,:),'markeredgecolor',[0,0,0],'markersize',4)
            end
            plot3(plot_data(:,x_axis_index(jj)),plot_data(:,y_axis_index(jj)),ones(size(plot_data(:,y_axis_index(jj))))*5,'--','lineWidth',0.5,'color',[127,127,127]./255);
        end
        xlabel(axis_labels{x_axis_index(jj)},'FontSize',18,'interpreter','latex');
        ylabel(axis_labels{y_axis_index(jj)},'FontSize',18,'interpreter','latex');
        
    end
    if contains(plot_type,'colour')==1
        if strcmp(input_name,'matSys')~=1
            hcb=colorbar('Ticks',[0 1],'TickLabels',colourmap_limits,'ticklabelinterpreter','latex','fontsize',14);
            colormap(colour_map);
            set(get(hcb,'Title'),'string',colourbar_label,'interpreter','latex');
        else
            h=legend(strrep(matSysCategories,"_"," "));
            set(h,'FontSize',11,'interpreter','latex','location','northeast');
        end
        
    end
    
    set(gca,'FontSize',18,'ticklabelinterpreter','latex')
    box off
    
    set(fig,'renderer','painters');
    
    %     set(gcf,'units','normalized','position',[0.2,0.6,0.5,0.8])
    
    
    if save_flag == 1
        plot_filepath=['results/MOBO/additive/' defect_type '/plots/' plot_type '/'];
        plot_name=[ plot_names{x_axis_index(jj)} '_' plot_names{y_axis_index(jj)} '_' input_name];
        try
            savefig([plot_filepath,plot_name]);
            saveas(gcf,[plot_filepath plot_name '.emf'])
        catch
            mkdir(plot_filepath);
            savefig([plot_filepath,plot_name]);
            saveas(gcf,[plot_filepath plot_name '.emf'])
        end
    end
end



end

function matSysCategories=create_matSysSategories()

Fibre_types=["E-Glass","S-Glass","PAN-CF","Pitch-CF"];

matSysCategories=[];
for ii=1:4
    for jj=1:4
        matSysCategories=[matSysCategories; strjoin([Fibre_types(ii) ", " Fibre_types(jj)],"")];
    end
    matSysCategories=[matSysCategories; strjoin(["non-hybrid " Fibre_types(ii)],"")]; % add non-hybrid cases
end



matSysCategories=categories(categorical(cellstr(matSysCategories)));
end