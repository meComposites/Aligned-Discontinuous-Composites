function [Final_Xfailure, Final_critical_strain_increment, Final_critical_cluster, RVE_data] = SFCfracture_Clusters_EffectsOfRandomness_bwconncomp_othernewX(deltaEpsilon,Ec,Eg,Vc,nc,Gcrit_c,Gcrit_g,n_rows,n_columns,Xfibres,Df,dc,dg,Vf,Xmax,iXmax,Enom,RVE_data,tm,lf,Giicm_CC,Sm_CC,gm_CC,Giicm_GC,Sm_GC,gm_GC,Giicm_GG,Sm_GG,gm_GG)
% dbstop if error

toughness_method='average toughness method';

init_num_clusters=100;

% Fibre_stiffness=bsxfun(@rdivide,Xfibres(:,:,2),deltaEpsilon);
Fibre_stiffness=RVE_data.Pre_defect_stiffness;

RVE_fibre_type=RVE_data.RVE_fibre_type;

%% set the fibre effectiveness increment

Damage_increment=0.01;
% Max_fibre_effectiveness=0.005;
Max_damage_threshold=1;
Damage_threshold=[Max_damage_threshold:-Damage_increment:0]';
Num_damage_increments=Max_damage_threshold/Damage_increment;

dim_d=Num_damage_increments-1;
dim_i=size(Xfibres,1);
dim_j=size(Xfibres,2);
dim_s=size(RVE_data.Stress_strain_unprocessed,1);
dim_k=init_num_clusters;

Initial_fibre_indices=[1:dim_i*dim_j];

Perfect_fibre_stresses=bsxfun(@times,repmat(Fibre_stiffness,1,1,dim_s),repmat(reshape(([1:dim_s]-1)*deltaEpsilon,1,1,dim_s),dim_i,dim_j));
Damage_state=bsxfun(@rdivide,(Perfect_fibre_stresses-Xfibres(:,:,1:dim_s)),Perfect_fibre_stresses);
%%%% add in failed fibres here at start of analysis if looking at flaws %%%
Damage_state(:,:,1)=zeros(dim_i,dim_j);

% fibres that never carry any stress are considered fully failed
[i_prefailed,j_prefailed]=find(bsxfun(@eq,Xfibres(:,:,2),0));
if isempty(i_prefailed)==0
    for i=1:length(i_prefailed)
        Damage_state(i_prefailed(i),j_prefailed(i),:)=1;
    end
    dim_k=dim_k+length(i_prefailed);
end

% calculate total fibre area (used later)
Fibre_areas=(RVE_fibre_type==1).*(pi()*(dc/2)^2)+(RVE_fibre_type==0).*(pi()*(dg/2)^2);
Total_fibre_area=pi()*sum(sum(RVE_fibre_type*(dc/2)^2+(1-RVE_fibre_type)*(dg/2)^2));
% Applied_stress=RVE_data.Stress_strain_unprocessed(:,2);

% Xcrit=zeros(dim_s,dim_d,dim_k);

%% loop through strain increment
reserve_factor=inf(dim_s,dim_d,dim_k);
for s=2:dim_s
    % loop through fibre effectiveness levels
    
    Xfibres_at_s=Xfibres(:,:,s);
    Damage_state_at_s=Damage_state(:,:,s);
    
    if s==40
        1;
    end
    
    ClustersAtDamageThreshold=zeros(dim_i,dim_j,dim_d);
    n_k=zeros(dim_d,1);
%     Xcrit=zeros(dim_d,dim_k);
%     Xcrit_record=zeros(dim_d,dim_k);
    
    % determine which fibres are within the cluster
%     FibresWithinCluster=bsxfun(@ge,repmat(reshape(Damage_state(:,:,s),dim_i,dim_j),1,1,dim_d),repmat(reshape(Damage_threshold(1:end-1),1,1,dim_d),dim_i,dim_j,1));
    FibresWithinCluster=bsxfun(@ge,repmat(reshape(Damage_state(:,:,s),dim_i,dim_j),1,1,dim_d),repmat(reshape(Damage_threshold(1:end-2),1,1,dim_d),dim_i,dim_j,1));
    
    % loop through damage threshold (looped this way to speed up)
    for d = 1:(dim_d)
        
%         [ClustersAtDamageThreshold(:,:,d) n_k(d)]=bwlabel(FibresWithinCluster(:,:,d),4);

        Cluster_conn_data=bwconncomp(FibresWithinCluster(:,:,d),4);
        
%         if n_k(d)>0
        if Cluster_conn_data.NumObjects>0
            
            while Cluster_conn_data.NumObjects>dim_k
                dim_k=dim_k+100;
                reserve_factor=cat(3,reserve_factor,inf(dim_s,dim_d,dim_k-size(reserve_factor,3)));
%                 Xcrit=cat(3,Xcrit,zeros(dim_s,dim_d,dim_k-size(Xcrit,3)));
%                 Xcrit_record=cat(2,Xcrit_record,zeros(dim_d,dim_k-size(Xcrit_record,2)));
            end
            
            Cluster_indices_record=[];

            for k=1:Cluster_conn_data.NumObjects
                %% fracture toughness calcs

%                 Cluster_indices=find(bsxfun(@eq,ClustersAtDamageThreshold(:,:,d),k));
                Cluster_indices=Cluster_conn_data.PixelIdxList{k};
                
%                 Residual_stress=sum(sum(bsxfun(@times,Xfibres_at_s(Cluster_indices),Fibre_areas(Cluster_indices))))/sum(sum(Fibre_areas(Cluster_indices)));
                
                %                 Residual_stress=sum(sum(bsxfun(@times,clusters,bsxfun(@times,Xfibres(:,:,s),Fibre_areas))))/sum(sum(bsxfun(@times,clusters,Fibre_areas)));
%                 Xcrit(d,k)=Applied_stress(s)-Residual_stress;
                
%                 Xcrit_record(d,k)=Xcrit(d,k);
                
                % only tensile loads will cause a crack to open, so ignore
                % compressive loads
                
                %                 Xcrit(Xcrit<0)=0;
                
                % only tensile loads will cause a crack to open, so ignore
                % compressive loads
%                 if Xcrit(d,k)<0
%                     Xcrit(d,k)=0;
%                     Cluster_indices=[];
%                 end
                
                Cluster_indices_record=[Cluster_indices_record; Cluster_indices];
                                
            end
            
            Applied_stress_fibre_indices=setdiff(Initial_fibre_indices,Cluster_indices_record);
            
            Applied_stress_fibre_area=sum(sum(Fibre_areas(Applied_stress_fibre_indices)));
            
%             Applied_stress=sum(sum(bsxfun(@times,Xfibres_at_s(Applied_stress_fibre_indices),Fibre_areas(Applied_stress_fibre_indices))))./Applied_stress_fibre_area
%             Applied_stress=sum(sum(bsxfun(@times,Xfibres_at_s(Applied_stress_fibre_indices),Fibre_areas(Applied_stress_fibre_indices))))./Total_fibre_area;
            
%             Residual_stress=sum(sum(bsxfun(@times,Xfibres_at_s(Cluster_indices_record),Fibre_areas(Cluster_indices_record))))./sum(sum(Fibre_areas(Cluster_indices_record)));
%             Residual_stress=sum(sum(bsxfun(@times,Xfibres_at_s(Cluster_indices_record),Fibre_areas(Cluster_indices_record))))./Total_fibre_area;
            
            Applied_force=sum(sum(bsxfun(@times,Xfibres_at_s(Applied_stress_fibre_indices),Fibre_areas(Applied_stress_fibre_indices))));
            Residual_force=sum(sum(bsxfun(@times,Xfibres_at_s(Cluster_indices_record),Fibre_areas(Cluster_indices_record))));
                       
            
%             Xcrit=max([Applied_stress-Residual_stress,0]);
            Xcrit=max([(Applied_force-Residual_force)./Total_fibre_area,0]);
%             Xcrit=max([(Applied_force-Residual_force)./Applied_stress_fibre_area,0]);
            
%             Broken_fibre_area=sum(sum(Fibre_areas(Cluster_indices_with_positive_driving_stress)));
%             Broken_fibre_area=sum(sum(bsxfun(@times,Damage_state_at_s(Cluster_indices_record),Fibre_areas(Cluster_indices_record))));
            
%             fncrit=(Total_fibre_area-Broken_fibre_area)/Total_fibre_area;
            
            for k=1:Cluster_conn_data.NumObjects
                
%                 Cluster_indices=find(bsxfun(@eq,ClustersAtDamageThreshold(:,:,d),k));
                Cluster_indices=Cluster_conn_data.PixelIdxList{k};
                
%                 Residual_stress=sum(sum(bsxfun(@times,Xfibres_at_s(Cluster_indices),Fibre_areas(Cluster_indices))))./Applied_stress_fibre_area;
%                 Residual_stress=sum(sum(bsxfun(@times,Xfibres_at_s(Cluster_indices),Fibre_areas(Cluster_indices))))./sum(sum(Fibre_areas(Cluster_indices)));
                
%                 Xcrit=max([Applied_stress-Residual_stress,0]);
                
%                 n_c_cluster=sum(sum(RVE_data.RVE_fibre_type(Cluster_indices)));
%                 n_g_cluster=sum(sum(1-RVE_data.RVE_fibre_type(Cluster_indices)));
                
                switch toughness_method
                    case 'surrounding fibres method'
                        clusters=bsxfun(@eq,ClustersAtDamageThreshold(:,:,d),k);
                        
                        % create a matrix of fibres that surround the cluster, by
                        % transforming the clusters matrix around itself (co-ordination
                        % number 4)
                        %                 surrounding_fibres_outer=[zeros(1,n_columns+2);zeros(n_rows,1),clusters,zeros(n_rows,1);zeros(1,n_columns+2)];
                        surrounding_fibres_outer=zeros(n_rows+2,n_columns+2);
                        surrounding_fibres_outer(2:end-1,2:end-1)=clusters;
                        surrounding_fibres_outer(1:end-2,2:end-1)=surrounding_fibres_outer(1:end-2,2:end-1)+clusters;
                        surrounding_fibres_outer(3:end,2:end-1)=surrounding_fibres_outer(3:end,2:end-1)+clusters;
                        surrounding_fibres_outer(2:end-1,1:end-2)=surrounding_fibres_outer(2:end-1,1:end-2)+clusters;
                        surrounding_fibres_outer(2:end-1,3:end)=surrounding_fibres_outer(2:end-1,3:end)+clusters;
                        surrounding_fibres_outer=(surrounding_fibres_outer>0);
                        surrounding_fibres=surrounding_fibres_outer(2:end-1,2:end-1);
                        
                        % create matrix of fibres which surround clusters within
                        % cross-section
                        surrounding_fibres=surrounding_fibres-clusters;
                        
                        nc_border=sum(sum(bsxfun(@times,surrounding_fibres,RVE_fibre_type)));
                        ng_border=sum(sum(bsxfun(@times,surrounding_fibres,(1-RVE_fibre_type))));
                        
                        Giicm=(nc_border*(dc^2)*Gcrit_c+ng_border*(dg^2)*Gcrit_g)/(nc_border*(dc^2)+ng_border*(dg^2));
                        
                    case 'average toughness method'
                        Giicm=Vc*Gcrit_c+(1-Vc)*Gcrit_g;
                end
                
                %% strain energy release rate calcs
                % FJ - aeq ~ a_crit (Eqn 17 in {S.Pimenta}) with some scaling?
                
                % equivalent crack radius: find equivalent radius of a circular
                % area containing fibre and matrix (modified from Eqn 17 in {S.Pimenta})
%                 aeq=0.5*sqrt((n_c_cluster*(dc^2)+n_g_cluster*(dg^2))/Vf);
                aeq=0.5*sqrt((sum(sum(bsxfun(@times,Damage_state_at_s(Cluster_indices),bsxfun(@times,RVE_data.RVE_fibre_type(Cluster_indices),dc^2))))+sum(sum(bsxfun(@times,Damage_state_at_s(Cluster_indices),bsxfun(@times,(1-RVE_data.RVE_fibre_type(Cluster_indices)),dg^2)))))/Vf);
                
                % Calculate J-integral
%                 Jpp=real(Xmax^2/Enom * 32/pi^3 * aeq .* (log( 1./cosd(90.* Xcrit(d,k)./(Xmax*fncrit) ) )));
                Jpp=real(Xmax^2/Enom * 32/pi^3 * aeq .* (log( 1./cosd(90.* Xcrit./(Xmax) ) )));
                
                % Calculate reserve factor
                reserve_factor(s,d,k)=Giicm/Jpp;
                
            end
        end
    end
    
    
    %     reserve_factor(reserve_factor==0)=Inf;
    
    if min(min(reserve_factor(s,:,:)))<1
        [reserve_factor_temp d_crit]=min(reserve_factor(s,:,:));
        [min_reserve_factor critical_cluster_index]=min(reserve_factor_temp);
        d_crit=d_crit(critical_cluster_index);
        RVE_data.Crit_damage_threshold=Damage_threshold(d_crit);
        Final_critical_strain_increment=s;
        Final_critical_cluster=bwlabel(FibresWithinCluster(:,:,d_crit),4)==critical_cluster_index;
        RVE_data.Cluster_size_record=sum(sum(Final_critical_cluster));
        Final_Xfailure = RVE_data.Stress_strain_unprocessed(s,2);
        RVE_data.Fracture_condition='True';
        RVE_data.Damage_state_at_failure=Damage_state(:,:,s);
        break
    end
end

if min(min(min(reserve_factor)))>=1
    RVE_data.Fracture_condition='False';
    [Final_Xfailure Final_critical_strain_increment] = max(RVE_data.Stress_strain_unprocessed(:,2));
    Final_critical_cluster=[];
    RVE_data.Crit_damage_threshold=[];
    RVE_data.Cluster_size_record=n_rows*n_columns;
    RVE_data.Damage_state_at_failure=Damage_state(:,:,Final_critical_strain_increment);
end

% disabled for FAST version
% store results in RVE data structure
% RVE_data.Reserve_factor_record=reserve_factor(:,:,Crit_fibre_effectiveness_index);

% RVE_data.Clusters_record=ClustersAtStrain(:,:,:,crit_fibre_effectiveness_index);
% RVE_data.Cluster_size_record=cluster_size(:,:,Crit_fibre_effectiveness_index);
% RVE_data.Min_crit_strain_all_fibre_effectiveness=critical_strain_increment*deltaEpsilon;
% RVE_data.CritStrainClustersAtEffectiveness=CritStrainClustersAtEffectiveness;
% RVE_data.Crit_fibre_effectiveness_index=Crit_fibre_effectiveness_index;
% RVE_data.CritClusterAtEffectiveness=CritClusterAtEffectiveness;
end