function [T1v,T2v,T1h,T2h] = aux_mapping_Tf_to_TI(Tf)
% Mapping fibre thicknesses to interactions

%vertical I, fibre above interaction, neighbour 3:
T1v=Tf(1:end-1,:,3);

%vertical I, fibre below interaction, neighbour 1:
T2v=Tf(2:end,:,1);

%horizontal I, fibre on the left interaction, neighbour 2:
T1h=Tf(:,1:end-1,2);

%horizontal I, fibre on the right interaction, neighbour 4:
T2h=Tf(:,2:end,4);

end

