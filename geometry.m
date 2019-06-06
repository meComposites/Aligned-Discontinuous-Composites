function [l,TA,TB,tm] = geometry(AT)

%--------------------------------------------------------------------------
% Geometry shared by Analytical model and FEA generator -------------------
%--------------------------------------------------------------------------

% Global geometric parameters:

l=20;        % Characteristic length (mm) l = 4.L
ET=0.2;     % Characteristic thickness (mm)
tm=0.01;    % Matrix thickness (mm)

l=2;        % Characteristic length (mm) l = 4.L
ET=0.007;     % Characteristic thickness (mm)
tm=0.001;    % Matrix thickness (mm)

%l=0.2;       % Spider silk - Bamboo
%ET=0.02;     % Spider silk - Bamboo
%tm=0.001;    % Spider silk - Bamboo

% Hybridisation (Thickness TA and TB):

if nargin == 0
    AT=0;   % Default hybridisation (same thickness)
end

TA=(1-AT)*ET/2;
TB=(1+AT)*ET/2;

end