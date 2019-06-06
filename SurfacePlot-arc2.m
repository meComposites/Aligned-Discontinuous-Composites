function [] = SurfacePlot( Lall,event )
% The function aims at representing the evolution of shear stress with both
% shear strength and position along the platelet. To do so, functions have
% been modified to collect the list of gCrit, and a list of point is
% created from the matrices field.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Recuperation of results via SL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(Lall)>1
    msgbox('SurfacePlot cannot be used for more than one input length')
elseif event<=0
    msgbox('event must be a positif integer')
else
    [~,~,~,~,fieldx,~,fieldS,~,~,gCrit]=SL(Lall);
    [~,~,~,~,~,~,~,Ng] = Inputs();
    
    NbrGCrit=size(gCrit,1);
    
    VecG=[];    % Creation of a matrix containing the steps of gamma
    for i=1:NbrGCrit-1
        VecG=[VecG linspace2(gCrit(i), gCrit(i+1), Ng+1)];
    end
    
                                        % Creation fo a matrix containing
    VecG=repmat(VecG,size(fieldx,2),1); % the value of gamma used to
    VecG=100*VecG';                     % calculate tau (gamma critical +
                                        % interpolation with Ng values)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Interpolation of results with scatteredInterpolant %%%%%%%%%%%%%%%%%%

    %%%%% creation of on x,y %%%%%%%%%%%%%%%
    
    if nargin == 1
        maxG=max(max(VecG));        % Maximum value for gamma
    else
        maxG=100*gCrit(min(size(gCrit,1),floor(event)+1));        % Maximum value for gamma zoomed
    end
    
    maxX=Lall;                  % Maximum value for x
    precisionG=maxG/1000;       % Interpolation step for gamma
    precisionX=maxX/100;        % Interpolation step for x
        
    xi=-maxX:precisionX:maxX;   % Creation of x scale
    yi=0:precisionG:maxG;       % Creation of gamma scale
    
    [x,y] = meshgrid(xi,yi);    % Creation of the grid
    
    %%%%% Creation of a list of points (x,g,tau) %%%%%%%%%%%%%%%
    
    ListPoints=[];
    for j=1:size(fieldx,2)
        for i=1:size(fieldx,1)
            if ~isnan(fieldS(i,j))      % if the point is activated (NaN)
                ListPoints=[ListPoints ;...
                    [fieldx(i,j) VecG(i,j) fieldS(i,j)]];
            end
        end
    end
    
    %%%%% Interpolation on (x,y) grid %%%%%%%%%%%%%%%
    
    %size(ListPoints)
    
    %%%%% Deletion of repeted lines in ListPoint
    [~,listunique]=unique(ListPoints(:,1:2),'rows');
    ListPoints=(ListPoints(listunique,:));
    %size(ListPoints)
    
    X=ListPoints(:,1);
    Y=ListPoints(:,2);
    S=ListPoints(:,3);
    
%     size(list2)
%     size(list)
%     figure
%     bar(list2)
%     figure
%     bar(sort(list))
    
    F = scatteredInterpolant(X,Y,S);    % Interpolation of (X,Y,S)
    InterpS=F(x,y);                     % Application on the grid (x,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     figure
     surface(x,y,InterpS,'LineStyle','none');
     colorbar
     title('Shear stress (MPa)')
     xlabel('x (mm)') 
     ylabel('Shear strain (%)')
     axis([-maxX,maxX,0,maxG])
     view(2)
     caxis([0,max(max(fieldS))])
end
end