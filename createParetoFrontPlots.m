%% create plain old Pareto fronts (no jazzy colours etc.) - use this to compare defective Paretos, or estimated Pareto fronts, etc.

function createParetoFrontPlots(paretoNames,paretoLegend,saveflag,markerColours,z_offsets,varargin)

xLabels={'Initial stiffness (GPa)','Ultimate strain (\%)','Pseudo-ductile strain (\%)'};
yLabels={'Ultimate strength (MPa)','Ultimate strength (MPa)','Yield strength (MPa)'};
designCases={'ultStrength_initS','ultStrength_ultStrain','yieldStrength_pdStrain'};
yLims=[0 2000;0 2000;0 2000];
xLims=[0 600;0 5;0 5];

if nargin~=4
    markerColours=[0,121,212;113,174,46;255,164,31;138,0,0;147,75,201]./255;
end

if nargin~=5
    z_offsets=1:length(paretoNames);
end

for designCaseIndex=1:3
    
    fig=figure();
    hold on
    
    for ii=1:length(paretoNames)
        load([paretoNames{ii} '.mat'])
        plot_data=eval(['paretoOptimal_' designCases{designCaseIndex}]);
        if designCaseIndex==1
            plot3(plot_data(:,1)./1000,plot_data(:,2),ones(size(plot_data,1),1).*z_offsets(ii),'o','markerfacecolor',markerColours(ii,:),'markeredgecolor',[0,0,0],'markersize',6)
        else
            plot3(plot_data(:,1),plot_data(:,2),ones(size(plot_data,1),1).*z_offsets(ii),'o','markerfacecolor',markerColours(ii,:),'markeredgecolor',[0,0,0],'markersize',6)
        end
        
    end
    
    xlabel(xLabels{designCaseIndex},'interpreter','latex','fontsize',18);
    xlim(xLims(designCaseIndex,:))
    ylabel(yLabels{designCaseIndex},'interpreter','latex','fontsize',18);
    ylim(yLims(designCaseIndex,:))
    set(gca,'ticklabelinterpreter','latex','fontsize',18);
    
    for ii=1:length(paretoNames)
        load([paretoNames{ii} '.mat'])
        plot_data=eval(['paretoOptimal_' designCases{designCaseIndex}]);
        if designCaseIndex==1
            plot3(plot_data(:,1)./1000,plot_data(:,2),ones(size(plot_data,1),1).*10,'--','lineWidth',0.5,'color',[127,127,127]./255);
        else
            plot3(plot_data(:,1),plot_data(:,2),ones(size(plot_data,1),1).*10,'--','lineWidth',0.5,'color',[127,127,127]./255);
        end
        
    end
    
    h=legend(paretoLegend);
    set(h,'interpreter','latex','fontsize',14,'location','northeast');
    
    set(fig,'renderer','painters'); % makes sure it's all vectorised
    
    switch saveflag
        case 'save'
            savefig(['results-second_pass/sensitivityAnalysis/sensitivity_Pareto_fronts/' paretoNames{ii} '_' designCases{designCaseIndex} '_Plot'])
            saveas(fig,[paretoNames{ii} '_' designCases{designCaseIndex} '_Plot'],'emf')
    end
end

end