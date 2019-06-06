function [T,tm,Eb,Nx,Ng] = SPInputs_VariableGiicSm(Ei,Ti,tmi)
% Defines Inputs for Shear-Lag model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry (characteristic length (L) is a variable for SL function)
% T=0.100;        % characteristic thickness, T (Figure 1.b)
T=Ti;        % characteristic thickness, T (Figure 1.b)
% tm=0.010;       % matrix thickness, t^m (Figure 1.b)
tm=tmi;       % matrix thickness, t^m (Figure 1.b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical properties (Table 5)
% Eb=100000;      % Stiffness modulus of the platelets, E^b
Eb=Ei;      % Stiffness modulus of the platelets, E^b
% FJ - Giic set to variable factor depending on fibre type (in
% librarygeneration function)
% Giicm=1;      % matrix's mode-II fracture toughness, G^m_IIc
% Giicm=2;      % matrix's mode-II fracture toughness, G^m_IIc
% Giicm=Giicm;      % matrix's mode-II fracture toughness, G^m_IIc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix constitutive law in shear 
% (from 0 to the maximum shear stress, Figure 2a)
% gmIn is g^[i], set of transition matrix shear strains
% SmIn is S^[i], set of transition matrix shear stresses 

% gmIn=[0.025,0.5,1]';      %nominal strain hardening:
% SmIn=[25,25,50]';         %nominal strain hardening:

%gmIn=[0.05]';            %nominal bilinear:
%SmIn=[50]';              %nominal bilinear:

%gmIn=[0.05,2.00]';       %nominal perfectly plastic, tm=0.010:
%SmIn=[50,50]';           %nominal perfectly plastic, tm=0.010:

%gmIn=[0.05,1.015]';      %nominal perfectly plastic, tm=0.020:
%SmIn=[50,50]';           %nominal perfectly plastic, tm=0.020:

% FJ - matrix strength allowed to vary between fibre types
% gmIn=[100/1500]';      %nominal strain hardening:
% SmIn=[100]';         %nominal strain hardening:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical variables
Nx=10;           % number of x-points within a subdomain for calculating 
                % stress and strain fields (x in Figure 2.b)
                
Ng=10;           % number of g0-points in between transition strains
                % (Ng+1 points for each set s in Figure 3)

end