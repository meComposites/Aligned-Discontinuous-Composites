function ductilityStudy_propertyMaps(batch_index)
load('pristine_VTF_data.mat')
index=find(Pseudo_ductile_strain==max(Pseudo_ductile_strain));

if batch_index>=1 && batch_index<=100
        % bdiamater only study (two HMC fibres)
        dc=[3:12]*0.001;
        dg=[3:12]*0.001;
        
        ii=ceil(batch_index/length(dg));
        jj=mod(batch_index,length(dg))+(mod(batch_index,length(dg))==0)*length(dg);
        
        Specimen_dataset = Specimen_Microscale_General(['diameterMap_dc_' num2str(dc(ii)) '_dg_' num2str(dg(jj))],Vf(index),Glass_fibre(index),EgAvg(index),EgCoV(index),Sg(index),mgWei(index),lgWei(index),dc(ii),Glass_fibre(index),EgAvg(index),EgCoV(index),Sg(index),mgWei(index),lgWei(index),dg(jj),1-Vc(index),lf(index)*1000,SmAvg(index),G(index),GiicmAvg(index),0.15,10,50000,80,80,0,0,0,0,0,0,0,1,0,1);
        
elseif batch_index>=101 && batch_index<=200
        % modulus only study (two HMC fibres)
        batch_index=batch_index-100;
    
        EcAvg=[73:(940-73)/9:940]*1000;
        EgAvg=[73:(940-73)/9:940]*1000;
        
        ii=ceil(batch_index/length(EgAvg));
        jj=mod(batch_index,length(EgAvg))+(mod(batch_index,length(EgAvg))==0)*length(EgAvg);
        
        Specimen_dataset = Specimen_Microscale_General(['modulusMap_EcAvg_' num2str(EcAvg(ii)) '_EgAvg_' num2str(EgAvg(jj))],Vf(index),Glass_fibre(index),EcAvg(ii),EgCoV(index),Sg(index),mgWei(index),lgWei(index),dg(index),Carbon_fibre(index),EgAvg(jj),EcCoV(index),Sc(index),mcWei(index),lcWei(index),dc(index),1-Vc(index),lf(index)*1000,SmAvg(index),G(index),GiicmAvg(index),0.15,10,50000,80,80,0,0,0,0,0,0,0,1,0,1);
        
elseif batch_index>=201 && batch_index<=300
        
        batch_index=batch_index-200;
    
        dc=[3:12]*0.001;
        EcAvg=[73:(940-73)/9:940]*1000;
        
        ii=ceil(batch_index/length(EcAvg));
        jj=mod(batch_index,length(EcAvg))+(mod(batch_index,length(EcAvg))==0)*length(EcAvg);
        
        Specimen_dataset = Specimen_Microscale_General(['structuralStiffnessMap_dg_' num2str(dc(ii)) '_EgAvg_' num2str(EcAvg(jj))],Vf(index),Glass_fibre(index),EgAvg(index),EgCoV(index),Sg(index),mgWei(index),lgWei(index),dg(index),Carbon_fibre(index),EcAvg(jj),EcCoV(index),Sc(index),mcWei(index),lcWei(index),dc(ii),1-Vc(index),lf(index)*1000,SmAvg(index),G(index),GiicmAvg(index),0.15,10,50000,80,80,0,0,0,0,0,0,0,1,0,1);
        
else
    error('incorrect batch_index')
end
end