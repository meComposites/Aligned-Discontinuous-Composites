load('D:\Box Sync\HiPerDuCT PhD\Modelling\Effect of randomness - personal folder\HPC_files\HPC_output_folder - important data from external drive\Defects\Individual influence of each defect\old_plots\defect_data_record.mat')
markerColors=[1,164/255,31/255;113/255,174/255,46/255;113/255,174/255,46/255];
markerStyle={'x','s','d'};
lineStyle={'-.','--','-'};
yLabels={'\% change in ultimate strain','\% change in ultimate strength','\% change initial stiffness'};
propertyNames={'ultStrain','ultStrength','initS'};

for property=[1,3,5]
    figure()
    hold on
    for matSys=[1,3]
        baseline=mean(mean(mean(data_record(1,:,matSys,property,:))));
        
        for defectType=1:3
            %             defectPercent=repmat([1:6]-1*5,1,3);
            defectPercent=repmat(([1:6]-1)*5,1,1);
            yAxis=mean(data_record(:,:,matSys,property,defectType)./baseline,2)*100;
            
            %                 plot((defectPercent-1)*5,mean(mean(data_record(defectPercent,:,matSys,property,defectType)./baseline)),'x','marker',markerStyle{defectType},'linewidth',3,'markersize',10,'color',markerColors(matSys,:));
            plot(defectPercent,yAxis,'-','linestyle',lineStyle{defectType},'linewidth',2,'color',markerColors(matSys,:));
            
        end
        
    end
ylabel(yLabels{(property-1)/2+1},'FontSize',14,'interpreter','latex');
xlabel('\% defect','FontSize',14,'interpreter','latex');
set(gca,'FontSize',14,'ticklabelinterpreter','latex')
box off
h=legend('6mm HSC w/ fragmentation defects','6mm HSC w/ interface defects','6mm HSC w/ fibre vacancy defects','6mm HMC/EG w/ fragmentation defects','6mm HMC/EG w/ interface defects','6mm HMC/EG w/ fibre vacancy defects');
set(h,'FontSize',12,'interpreter','latex','location','northeast');
xlim([0 25]);    
ylim([40 150]);
    
saveas(gcf,[propertyNames{(property-1)/2+1} '_percentage_defects.emf'])
savefig([propertyNames{(property-1)/2+1} '_percentage_defects'])

end