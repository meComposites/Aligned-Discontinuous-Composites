function [av,bv,ah,bh] = aux_mapping_f_to_I(f)
% Mapping fibre discontinuities to interactions

%vertical I, fibre above interaction, neighbour 3:
av=f(1:end-1,:,:);

%vertical I, fibre below interaction, neighbour 1:
bv=f(2:end,:,:);

%horizontal I, fibre on the left interaction, neighbour 2:
ah=f(:,1:end-1,:);

%horizontal I, fibre on the right interaction, neighbour 4:
bh=f(:,2:end,:);

end

