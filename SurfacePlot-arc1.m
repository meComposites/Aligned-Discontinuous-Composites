function [] = SurfacePlot( Lall )
% The function aims at representing the evolution of shear stress with both
% shear strength and position along the platelet. To do so, functions have
% been modified to collect the list of gCrit, and a list of point is
% created from the matrices field.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Recuperation of results via SL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(Lall)>1
    msgbox('SurfacePlot cannot be used for more than one input length')
else
    [~,~,~,~,fieldx,~,fieldS,~,~,gCrit]=SL(Lall);
    [~,~,~,~,~,~,~,Ng] = Inputs();
    
    NbrGCrit=size(gCrit,1);
    
    VecG=[];    % Creation of a matrix containing the steps of gamma
    for i=1:NbrGCrit-1
        VecG=[VecG linspace2(gCrit(i), gCrit(i+1), Ng+1)];
    end
    
    %VecG=unique(VecG);                  % Creation fo a matrix containing
    %VecG=VecG(1,1:size(VecG,2)-1);      % the value of gamma used to
    VecG=repmat(VecG,size(fieldx,2),1); % calculate tau (gamma critical +
    VecG=100*VecG';                     % interpolation with Ng values)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Interpolation of results with scatteredInterpolant %%%%%%%%%%%%%%%%%%

    %%%%% creation of on x,y %%%%%%%%%%%%%%%
    
    precisionG=0.0001;             % Interpolation step for gamma
    precisionX=0.01;            % Interpolation step for x
    maxG=max(max(VecG));        % Maximum value for gamma
    maxG=1;        % Maximum value for gamma
    maxX=Lall;                  % Maximum value for x
    
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
    
    %size(ListPoints,1)
    
    %%%%% Interpolation on (x,y) grid %%%%%%%%%%%%%%%
    
    X=ListPoints(:,1);
    Y=ListPoints(:,2);
    S=ListPoints(:,3);
    
    figure
    scatter3(X,Y,S,30,S,'fill')
    view(2)
    
%     writerObj = VideoWriter('test.avi');
%     writerObj.FrameRate = 3;
%     open(writerObj);
%     for k = 1:50
% 	axis([-4 4 -0.2*2^(min(25,k))*10^(-5) 1.5*2^(min(25,k))*10^(-5)])
% % 	M(k) = getframe;
%     frame=getframe;
%     writeVideo(writerObj,frame);
%     end
    
%     close(writerObj);
    
%     figure
%     movie(M,1,5)
%     ListPoints;
    
    %G=TriScatteredInterp(X,Y,S);
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
     
         writerObj = VideoWriter('test2.avi');
    writerObj.FrameRate = 3;
    open(writerObj);
    for k = 1:50
	axis([-4 4 -0.2*2^(min(25,k))*10^(-5) 1.5*2^(min(25,k))*10^(-5)])
% 	M(k) = getframe;
    frame=getframe;
    writeVideo(writerObj,frame);
    end
    close(writerObj);
    
end
end