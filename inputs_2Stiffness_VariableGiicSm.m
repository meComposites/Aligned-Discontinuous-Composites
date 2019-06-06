function [G,Sm,gm,Y] = inputs_2Stiffness_VariableGiicSm(EE,DE,ET,DT,tm,Gic,MLawS,MLawg)

%--------------------------------------------------------------------------
% Inputs ------------------------------------------------------------------
%--------------------------------------------------------------------------


% Matrix constitutive law in shear: ---------------------------------------

% Strain Hardening %
% MLawS=[40       40      80];        % Interface shear stress
% MLawg=[40/1500  0.5     1];   % Interface shear strain

% Trilinear %
% MLawS=[40         80];        % Interface shear stress
% MLawg=[40/1500    .5];   % Interface shear strain

% % Bilinear %
% MLawS=[80];        % Interface shear stress
% MLawg=[80/1500];   % Interface shear strain
% Gic=2;                   % Toughness
% 
% % Experimental guess
% MLawS=[50];        % Interface shear stress
% MLawg=[50/2700];   % Interface shear strain
% Gic=1;                   % Toughness

% Experimental guess
% MLawS=[100];        % Interface shear stress
% MLawg=[100/1500];   % Interface shear strain
% Gic=1;             % Toughness
% Gic=2;             % Toughness

MLawG=MLawS(1)/MLawg(1);
for im=1:length(MLawg)-1
    MLawG=[MLawG (MLawS(im+1)-MLawS(im))/(MLawg(im+1)-MLawg(im))];
end
[Sm,gm,G]=findLastg(MLawS,MLawg,MLawG,Gic,tm);

% Pre-calc ----------------------------------------------------------------

for ig=1:1:length(G)
    Y(ig)=sqrt(8*abs(G(ig))/tm*(EE*ET+DE*DT)/((EE^2-DE^2)*(ET^2-DT^2)));
end

end