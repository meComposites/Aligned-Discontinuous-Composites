function Gcrit = findGcrit(lf,Vf,Gm,Tu0,dc,alphac,Sm,lcWei,XcWei,mcWei,paramStudy,plotFig)


% --- Definition of spatial step for numerical integration ---
% dl=0.01;
dl=lf/300;

if paramStudy==1
    Gic=[0];
    lic=[0];
    Nic=[0];
    
    dlic=0.01;
    for li=dlic:dlic:lf
        Gcritl=findGcritSingleLength(li,Vf,Gm,Tu0,dc,alphac,Sm,lcWei,XcWei,mcWei,paramStudy);
        Gic=[Gic Gcritl];
        lic=[lic li];
    end
    
    % --- Plot Gic(lic)
    if plotFig==1
        figure
        plot(lic,Gic)
    end
    
    Gcrit=Gic(end);

else
    Gcrit=findGcritSingleLength(lf,Vf,Gm,Tu0,dc,alphac,Sm,lcWei,XcWei,mcWei,paramStudy);
end

1;
    function Gcritl = findGcritSingleLength(li,Vf,Gm,Tu0,dc,alphac,Sm,lcWei,XcWei,mcWei,paramStudy)
        
        % precalculation
        ldeb=dl:dl:li/2;
        Tu=Tu0*(1-Vf);
        gdeb=2*Vf*Gm/dc*ldeb;
        gpo=Vf*Tu/dc*ldeb.^2;
        lo=ones(1,length(ldeb))*1/(li/2);
        %SLo=exp(-li./lcWei./alphac*((4*max(Sm)*ldeb)/(dc*XcWei/gamma(1+1/mcWei))).^mcWei);
        SLo=exp(-ldeb./lcWei./alphac.*((4*max(Sm)*ldeb)/(dc*XcWei/gamma(1+1/mcWei))).^mcWei);
        
        % --- Definition of number of cracks ---
        Ncracks=4;
        %Ncracks=round(log2(lf/ldeb(sum(SLo>0.5))));
        %Ncracks=log2(10*lf);
        
        % --- Initialisation
        pdf{1}=lo;
        Mgdeb(1)=sum(gdeb.*(lo).*SLo.*dl);
        Mgpo(1)=sum(gpo.*(lo).*SLo.*dl);
        NewCrack{1}=lo'; % Survival probability of new (broken) overlaps only
        wba{1}=zeros(1,length(ldeb)); % Survival probability of new (broken) overlaps only
        % wba{1}=lo.*SLo;
        
        % --- Plot of initialisation (only for single length)
        if paramStudy==0
            if plotFig==1
                figure
                xlabel('length (mm)','interpreter','latex')
                ylabel('overlap distribution','interpreter','latex')
                set(gca,'TickLabelInterpreter', 'latex');
                set(gcf, 'PaperSize', [6.5 5.5]);
                set(gcf, 'Position', [2500 500 270 200])
                box off
                hold on
                plot(ldeb,lo)
                plot(ldeb,lo.*SLo)
            end
        end
        
        % --- Loop on number of cracks considered
        for ic=1:Ncracks
            %     err=1;
            %     ic=1;
            %     while abs(err)>0.01
            a=3/4;%3/4;
            
            % Survival probability of all overlaps
%             NewCrack{ic}=sum(bsxfun(@times,bsxfun(@times,bsxfun(@times,triu(ones(length(ldeb))),(1-SLo)),pdf{ic}),(1./ldeb)).*dl,2);
%             pdf{ic+1}=pdf{ic}.*SLo+a*NewCrack{ic}';
            
            % Survival probability of new (broken) overlaps only
            NewCrack{ic+1}=a*sum(bsxfun(@times,bsxfun(@times,bsxfun(@times,triu(ones(length(ldeb))),(1-SLo)),NewCrack{ic}'),(1./ldeb)).*dl,2);
            wba{ic+1}=wba{ic}+NewCrack{ic}'.*SLo;
            pdf{ic+1}=wba{ic+1}+NewCrack{ic+1}';
            
            % --- Plot for ic fibre breaks considered (only for single length)
            if paramStudy==0
                if plotFig==1
                    hold on
                    plot(ldeb,pdf{ic+1},'k--')
                end
            end
            
            Mgdeb(ic+1)=sum(gdeb.*(pdf{ic+1}).*dl);
            Mgpo(ic+1)=sum(gpo.*(pdf{ic+1}).*dl);
            
            %     err=(Mgdeb(ic+1)+Mgpo(ic+1)-(Mgdeb(ic)+Mgpo(ic)))/(Mgdeb(ic+1)+Mgpo(ic+1));
            %     ic=ic+1;
        end
        
        GcritlN=Mgdeb+Mgpo;
        Gcritl=GcritlN(end);
        
        % --- Plot for ic fibre breaks considered (only for single length)
        if paramStudy==0
            if plotFig==1
                figure
                bar([0:length(GcritlN)-1],GcritlN)
                xlabel('number of cracks','interpreter','latex')
                ylabel('$\mathcal{G}_\textrm{ic}$','interpreter','latex')
                % ylabel('$\frac{g_\textrm{deb}}{g_{\textrm{deb}_0}}$','interpreter','latex')
                set(gca,'TickLabelInterpreter', 'latex');
                set(gcf, 'PaperSize', [6.5 5.5]);
                set(gcf, 'Position', [2500 500 270 200])
                box off
            end
        end
        
    end




end