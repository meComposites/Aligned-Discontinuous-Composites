function [x,DX,S,g,nMax,colL,ijActive,gCrit] ...
    = SPSLResponse(N,type,Gm,y,gm,Sm,Lpz,L,T,tm,Eb,Nx,Ng)
% Main calculations for Shear-Lag Response

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initiating matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kmax=2*N;           % Calculates maximum number of transitions. 
if Lpz(N-1,1)>L     % There will be N transitions for activating subdomains
    kmax=kmax-1;    % 2:N+1, and N transitions for deactivating 
end                 % subdomains 1:N. If Lpz(N)>L, the activation of N+1 
                    % and deactivation of N will occur at the same moment.

nMax=sum(min(Lpz)<L)+1;         % Maximum n in the overlap 
                                % (maximum nMax+1 subdomains)

sMatrix=NaN(kmax+1,nMax+1);     % matrix with all sets of active subdomains
gCrit=NaN(kmax+1,1);            % vector of transition strains, \gamma^crit
x=NaN((Ng+1)*(kmax),(nMax+1)*(Nx+1));       % coordinate along overlap, x
S=x;                            % shear stress in the matrix, \gamma
g=x;                            % shear strain in the matrix, \tau
DX=x;                           % tensile stress difference, \Delta\sigma

colL=NaN((Ng+1)*(kmax),1);      % Vector with location of x=L fields
ijActive=NaN((Ng+1)*(kmax),2);  % matrix with subdomains i (centre) 
                                % and j=i+n (edge)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initiating for first domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gCrit(1)=0;
k=0;
sMatrix(1,1)=1;
ikplus1=1;      % subdomain in the centre for next set s_(k+1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop for main calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while(ikplus1~=N+1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k=k+1;
    sk=sMatrix(k,~isnan(sMatrix(k,:))); % new sk, without NaNs
    i=sk(1);
    n=length(sk)-1;

    % Calculating lines for increments of g0
    lines=(k-1)*(Ng+1)+1:k*(Ng+1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Select transition case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n==0         % if there is ONLY ONE SUBDOMAIN
        
        % Select transition case (only one subdomain)
        if i<N      % if it is NOT THE FINAL SUBDOMAIN  
            gCrit(k+1)=SPTransSing(type(i),Gm(i),y(i),...
                gm(i+1),Sm(i+1),L,T,tm,Eb);     % Activation of central 
            sMatrix(k+1,1:2)=[i,i+1];           % subdomain in Single Field
                    
        % Select transition case (only one subdomain)               
        else        % if it is THE FINAL SUBDOMAIN
            gCrit(k+1)=gm(end);                 % Final transition to
            sMatrix(k+1,1)=N+1;                 % fully debonded overlap
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate Fields for single subdomain (at the centre of overlap)
        cols=1:Nx+1;        % Update columns for fields
        
        li=L*ones(Ng+1,1);  % Calculate length of single fields (=L)
        
        % Calculate fields at the centre of the overlap
        DX0i=zeros(Ng+1,1);
        g0i=linspace(gCrit(k),gCrit(k+1),Ng+1)';
        S0i=interp1(gm,Sm,g0i,'linear');
        
        % Calculate full fields for single subdomain in the overlap
        [x(lines,cols),DX(lines,cols),S(lines,cols),g(lines,cols)]=...
        SPFieldsSing(type(i),Gm(i),y(i),g0i,S0i,DX0i,li,T,tm,Eb,Nx,Ng);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Select transition case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else            % if there are MULTIPLE subdomains

        % Select transition case (multiple subdomains)
        if sk(end)==N+1 || Lpz(sk(2)-1,n)>L     
            % If cannot add a new subdomain to the overlap
            gCrit(k+1)=gm(sk(1)+1);             % Deactivation of edge 
            sMatrix(k+1,1:n)=sk(2:end);         % subdomain (Table 4 Eq.a)
            
        % Select transition case (multiple subdomains)
        else% If can add next subdomain to the overlap
            gCrit(k+1)=SPTransMult(type,y,...     % Activation of 
                gm,Sm,sk,gCrit(k),L,T,tm,Eb);   % central subdomain in 
            sMatrix(k+1,1:n+2)=[sk,sk(end)+1];  % Multiple Field
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate Fields for multiple subdomains in the overlap %%%%%%%%%
        cols=1:(Nx+1)*(n+1);	% Update columns for fields
        
        % Calculate full fields for multiple subdomain in the overlap
        [x(lines,cols),DX(lines,cols),S(lines,cols),g(lines,cols)]=...
        SPFieldsMult(N,type,Gm,y,gm,Sm,sk,...
        gCrit(k),gCrit(k+1),L,T,tm,Eb,Nx,Ng);
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Storing general information about the current set of subdomains %%%%%
    ijActive(lines,1)=sk(1);     	% Central subdomain for current set;
    ijActive(lines,2)=sk(end);    	% Edge subdomain for next set;
    colL(lines,1)=(n+1)*(Nx+1);     % Column for x=L;
    ikplus1=sMatrix(k+1,1);       	% Central subdomain for next set.
    
end     % Return to main loop

ijActive(end+1,:)=N+1;      % Adds final step with i=N+1, i+n=N+1
end