%% creates plot of all sensitivity data to show which sensitivity cases most adversely affect the chosen Pareto front

function plot_all_Pareto_sensitivity_cases(defect_case,sensitivity_case,save_flag)
% load correct sensitivity data
dbstop if error
switch defect_case
    case 'pristine'
        switch sensitivity_case
            case 'no_small_hybrids'
                load('results-second_pass/estimatedParetoFronts/estimatedParetoFronts.mat'); 
                load('results-second_pass/sensitivityAnalysis/no_small_hybrids/paretoFrontSensitivityAnalysis.mat')
            case 'initial'
                load('results-second_pass/estimatedParetoFronts/estimatedParetoFronts.mat'); 
                load('results-second_pass/sensitivityAnalysis/initial_sensitivity_analysis/paretoFrontSensitivityAnalysis.mat')
        end
end

% sort Pareto front data
pareto_front_data=struct('paretoOptimal_ultStrength_initS',paretoOptimal_ultStrength_initS,'paretoOptimal_ultStrength_ultStrain',paretoOptimal_ultStrength_ultStrain,'paretoOptimal_yieldStrength_pdStrain',paretoOptimal_yieldStrength_pdStrain);
data_fields=fieldnames(pareto_front_data);

constraintViolations=reshape(mode(constraints<=0,2).*(-2)+1,size(meanProperties,1),11);
constraintViolations=mode(constraintViolations,2); % constraint violations for whole set

% get rid of baseline sensitivity case in centre of 3rd D of array
meanProperties=cat(3,meanProperties(:,:,1:5),meanProperties(:,:,7:end));
meanProperties=meanProperties(constraintViolations<0,:,:); 
MO_case=MO_case(constraintViolations<0);

% define which properties to use for each Pareto front
x_axis_index=[1,2,3];
y_axis_index=[4,4,5];
MO_cases=unique(MO_case);

% define limits and labels for plots
xlims=[0 600;0 5;0 5];
ylims=[0 2000];
xlabels={'Initial stiffness (GPa)','Ultimate strain (\%)','Pseudo-ductile strain (\%)'};
ylabels={'Ultimate strength (MPa)','Ultimate strength (MPa)','Yield strength (MPa)'};
sensitivity_labels={'$V_\mathrm{c}$','$l_\mathrm{f}$','$\textbf{E} \left(S_\mathrm{m}\right)$','$G$','$\mathcal{G}_\mathrm{II_m}$'};
sensitivity_markers={'v','^'};
marker_colours=colormap(lines);
marker_colours=marker_colours([1:5],:);


% create plots
for Pareto_index=1:3
    fig=figure();
    hold on
    for sensitivity_index=1:size(meanProperties,3)
        % select correct colour and marker style
        colour_index=mod(sensitivity_index,5);
        if colour_index==0
            colour_index=5;
        end
        marker_index=floor(sensitivity_index/6)+1;
%         if marker_index==0
%             marker_index=1;% make sure arrows pointing right direction
%         end
        
        % use correct plot data
        if Pareto_index==1
            plot_data=[meanProperties(strcmp(MO_case,MO_cases(Pareto_index)),x_axis_index(Pareto_index),sensitivity_index)./1000,meanProperties(strcmp(MO_case,MO_cases(Pareto_index)),y_axis_index(Pareto_index),sensitivity_index)];
        else
            plot_data=[meanProperties(strcmp(MO_case,MO_cases(Pareto_index)),x_axis_index(Pareto_index),sensitivity_index),meanProperties(strcmp(MO_case,MO_cases(Pareto_index)),y_axis_index(Pareto_index),sensitivity_index)];
        end
        
        % pot data, with random overlap to prevent all of same colour ontop
        plot3(plot_data(:,1),plot_data(:,2),rand(size(plot_data)),sensitivity_markers{marker_index},'markerfacecolor',marker_colours(colour_index,:),'markeredgecolor',[1,1,1],'markersize',8)
    end
    
    if Pareto_index==1
        plot3(pareto_front_data.(data_fields{Pareto_index})(:,1)./1000,pareto_front_data.(data_fields{Pareto_index})(:,2),ones(size(pareto_front_data.(data_fields{Pareto_index})))*5,'--','lineWidth',2,'color',[127,127,127]./255);
    else
        plot3(pareto_front_data.(data_fields{Pareto_index})(:,1),pareto_front_data.(data_fields{Pareto_index})(:,2),ones(size(pareto_front_data.(data_fields{Pareto_index})))*5,'--','lineWidth',2,'color',[127,127,127]./255);
    end
    
    %set limits etc
    xlim(xlims(Pareto_index,:));
    ylim(ylims);
    xlabel(xlabels{Pareto_index},'interpreter','latex','fontsize',18);
    ylabel(ylabels{Pareto_index},'interpreter','latex','fontsize',18);
    set(gca,'ticklabelinterpreter','latex','fontsize',18)
    box off
    
    % create snazzy legend
    nrow = 5;
    ncol = 2;
    ngrp = nrow*ncol;
    rowlabel = sensitivity_labels;
    collabel = {'$- 10\%$','$+ 10\%$'};
    [ypt, xpt] = ndgrid(1:nrow, 1:ncol);
    ypt(:,2)=ypt(:,2)+0.15;
    xpt(:,1)=xpt(:,1)+0.1725;
    xpt(:,2)=xpt(:,2)+0.35;
    pos = get(gca, 'position');
    sz = [0.28 0.27];
    legax = axes('position', [pos(1)+pos(3)-sz(1)-0.025 pos(2)+pos(4)-sz(2)-0.025 sz], ...
        'xlim', [-1 ncol+1], 'ylim', [-1 nrow+1], 'ydir', 'reverse', ...
        'xtick', [], 'ytick', [], 'box', 'on');
    hold on;
    for ii=1:2
        scatter(xpt(:,ii), ypt(:,ii),48,marker_colours,'filled',sensitivity_markers{ii},'markeredgecolor',[1 1 1]);
    end
    text([1.1725,2.35], zeros(1,ncol)-0.1, collabel, 'horiz', 'center','interpreter','latex','fontsize',12);
    text(ones(1,nrow)*0.5, [1:nrow], rowlabel, 'horiz', 'right','interpreter','latex','fontsize',12);
    
    % make sure still vectorised
    set(fig,'renderer','painters');
    
    % save figures
    switch save_flag
        case 'save'
            switch sensitivity_case
                case 'no_small_hybrids'
                    plot_filepath=['results-second_pass/sensitivityAnalysis/no_small_hybrids/allSensitivityCaseParetos/'];
                case 'initial'
                    plot_filepath=['results-second_pass/sensitivityAnalysis/initial_sensitivity_analysis/allSensitivityCaseParetos/'];
            end
%         plot_filepath='';
        plot_name=char(MO_cases(Pareto_index));
        
        savefig([plot_filepath,plot_name]);
        saveas(gcf,[plot_filepath plot_name '.emf'])
        
%         try
%             savefig([plot_filepath,plot_name]);
%             saveas(gcf,[plot_filepath plot_name '.emf'])
%         catch
%             mkdir(plot_filepath);
%             savefig([plot_filepath,plot_name]);
%             saveas(gcf,[plot_filepath plot_name '.emf'])
%         end
    end
    
end



end
