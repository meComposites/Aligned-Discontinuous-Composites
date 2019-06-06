function [segmentsI,xDiscI] = Interaction_segments(Lf,xDisca,xDiscb)
% Calculates matrix of segments of one interaction...
% based on discontinuities in fibres a and b

% Matrix of fibre discontinuities in the interaction
% 2nd row: 0=a, 1=b. Added 1st discontinuity outside RVE (+lf)
xDiscI=sortrows([   [xDisca(:)          zeros(length(xDisca),1)];...
                    [xDiscb(:)          ones(length(xDiscb),1)]]);
xDiscI=[xDiscI;     [xDiscI(1,1)+Lf     xDiscI(1,2)]];

% Matrix of segments in the interaction
% 2nd row: average of fibre indeces defining ends of segment
segmentsI=xDiscI(2:end,1)-xDiscI(1:end-1,1);
segmentsI(:,2)=(xDiscI(2:end,2)+xDiscI(1:end-1,2))./2;

end

