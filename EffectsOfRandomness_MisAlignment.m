function [Specimen_dataset] = EffectsOfRandomness_MisAlignment(Specimen_dataset,Constituent_properties,Alignment_distribution,Sigma,Alignment_method)
%% gather material data
% matrix
E11m=Constituent_properties.E11m;
E22m=Constituent_properties.E22m;
G12m=Constituent_properties.G12m;
v12m=Constituent_properties.v12m;

% glass fibre
E11g=Constituent_properties.E11g;
E22g=Constituent_properties.E22g;
G12g=Constituent_properties.G12g;
v12g=Constituent_properties.v12g;

% carbon fibre
E11c=Constituent_properties.E11c;
E22c=Constituent_properties.E22c;
G12c=Constituent_properties.G12c;
v12c=Constituent_properties.v12c;

% composite
Vf=Constituent_properties.Vf;
Vc=Constituent_properties.Vc;

%% Initial alignment distribution
switch Alignment_distribution
    case 'Yu2014'
        % initial angles from {Yu2014}
        pdfA=[-14 -12 -10 -8 -6 -4 -2 0 2 4 6 8 10 12];
        
        % initial angles distribution from {Yu2014}
        pdf=[0.007538152 0.005585776 0.005541936 0.048829431 0.050213025 0.132786514 0.172502589 0.264837957 0.20431619 0.063558419 0.02565602 0.007515476 0.005089178 0.0052823 ];
        
    case 'Sanadi1985'
        pdfA=-65:5:65;
        pdf=[0.00854199270560672,0.00979198159776840,0.0108900174800063,0.0124147855837546,0.0143929935092544,0.0166446560636039,0.0187454888879347,0.0208914188974543,0.0240076360646077,0.0295605553239190,0.0355001258743704,0.0426348746000149,0.162037997458971,0.182057209237711,0.175396688557475,0.0458793184709564,0.0327333241598710,0.0276430046284372,0.0236393281695494,0.0203873816179659,0.0179134212783292,0.0158136517712369,0.0138160080231605,0.0117955669124137,0.00987866665757605,0.00896321475306085,0.00802869171499021];
        
    case 'Normal_distribution'
        % normal distribution made to emulate alignment distribution
        x = [-65:1:65];
        norm = normpdf(x,0,Sigma);
        
        pdfA=x;
        
        pdf=norm/sum(norm);
        
end

%convert to rads
pdfA=pdfA*pi()/180;

%% perform initial alignment analysis to give misaligned stress-strain
% curve. The effects of fibre re-alignment can be calculated after
% collect the original stress-strain curve data
Orig_SS_curve=Specimen_dataset.Stress_strain_fracture;
Specimen_dataset.Misaligned_stress_strain(:,1)=Orig_SS_curve(:,1);
Specimen_dataset.Misaligned_stress_strain(1,2)=0;

% for every interval in the stress strain curve, evaluate the misaligned
% stiffness of the composite in the loading direction
for ii = 1:size(Orig_SS_curve,1)-1
    % collect the gradient of the stress-strain curve in this interval
    Orig_SS_grad(ii)=((Orig_SS_curve(ii+1,2)-Orig_SS_curve(ii,2))/(Orig_SS_curve(ii+1,1)-Orig_SS_curve(ii,1)))/10;
    
    % find the equivalent stiffness of the lamina, based on the rule of
    % mixtures
    E11=Orig_SS_grad(ii);
    E22=1/( (1-Vf)/E22m + (Vf*Vc)/E22c + (Vf*(1-Vc))/E22g );
    G12=1/( (1-Vf)/G12m + (Vf*Vc)/G12c + (Vf*(1-Vc))/G12g );
    v12=(1-Vf)*v12m + (Vf*Vc)*v12c + (Vf*(1-Vc))*v12g;
    
    E11=Orig_SS_grad(ii);
    E22=1/( (1-Vf)/E22m + (Vf*Vc)/E22c + (Vf*(1-Vc))/E22g );
    G12=1/( (1-Vf)/G12m + (Vf*Vc)/G12c + (Vf*(1-Vc))/G12g );
    v12=(1-Vf)*v12m + (Vf*Vc)*v12c + (Vf*(1-Vc))*v12g;
    
    % calculate reduced stiffness matrix terms
    Q11= E11^2/(E11-v12^2*E22);
    Q12= v12*E11*E22/(E11-v12^2*E22);
    Q22= E11*E22/(E11-v12^2*E22);
    Q66= G12;
    
    % loop through every entry in the alignment distribution
    for jj=1:length(pdfA)
        % convert to rads
        theta=pdfA(jj);
        
        % use these to simplify transformation matrix
        m=cos(theta);
        n=sin(theta);
        
        % create transformation matrix
        Mt=[m^4     n^4     2*m^2*n^2   4*m^2*n^2;
            n^4     m^4     2*m^2*n^2   4*m^2*n^2;
            m^2*n^2 m^2*n^2 m^4+n^4     -4*m^2*n^2];
        
        % find tansformed stiffness matrix
        Qt=Mt*[Q11 Q22 Q12 Q66]';
        
        % need to check what this does
        Qtm(1:3,jj)=Qt;
    end
    
    % find distributed tansformed stiffness values
    Q11b=sum(Qtm(1,:).*pdf);
    Q22b=sum(Qtm(2,:).*pdf);
    Q12b=sum(Qtm(3,:).*pdf);
    
    % calculate misaligned composite stiffness in x-direction
    Misaligned_SS_grad(ii)=(Q11b*Q22b-Q12b^2)/Q22b;
    
    % find misaligned SS curve
    Strain_increment=(Specimen_dataset.Misaligned_stress_strain(ii+1,1)-Specimen_dataset.Misaligned_stress_strain(ii,1));
    Specimen_dataset.Misaligned_stress_strain(ii+1,2)=Specimen_dataset.Misaligned_stress_strain(ii,2)+Misaligned_SS_grad(ii)*Strain_increment*10;
    
    % if realignment is included, change the angle of the fibres
    % accordingly
    switch Alignment_method
        case 'Realignment'
            pdfA=atan(tan(pdfA).*((1-v12*Strain_increment)/(1+Strain_increment)));       
    end
    
end



end