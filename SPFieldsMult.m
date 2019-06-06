function [x,DX,S,g]=SPFieldsMult(...
    N,type,Gm,y,gm,Sm,sk,g0min,g0max,L,T,tm,Eb,Nx,Ng)
% Calculates stress and strain fields for 
% several subdomains in the overlap length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initiating matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=length(sk)-1;

x=NaN(Ng+1,(n+1)*(Nx+1));
S=x; 
g=x;
DX=x;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fields for subdomain at the centre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=sk(1);                % Identifying subdomain at the centre

g0i=10.^linspace(log10(g0min),log10(g0max),Ng+1)';  % Calculating 
g0i(1)=g0min; g0i(end)=g0max;                       % fields
S0i=interp1(gm,Sm,g0i,'linear');                    % at
DX0i=zeros(Ng+1,1);                                 % x=0

col0=1;                 % Updating columns 
colL=Nx+1;              % for fields in the
cols=col0:colL;         % central subdomain

[li,~]=SPLengthSingCentre(...                         % Calculating length
    type(i),y(i),g0i,gm(i+1),S0i,Sm(i+1),T,tm,Eb);  % of subdomain i

% Calculating full fields for subdomain i at the centre of the overlap
[x(:,cols),DX(:,cols),S(:,cols),g(:,cols)]=...
    SPFieldsSing(type(i),Gm(i),y(i),g0i,S0i,DX0i,li,T,tm,Eb,Nx,Ng);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fields for completely developed subdomains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% (not at the centre nor at the edge)        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=sk(2:n);          % Identifying subdomain
 
    x0i=x(:,colL);                      % Calculating
    DX0i=DX(:,colL);                    % fields at       
    g0i=gm(i)*ones(Ng+1,1);             % x=x0(i)
    S0i=interp1(gm,Sm,g0i,'linear');    % (which corresponds to xL(i-1))

    col0=colL+1;        % Updating columns
    colL=col0+Nx;       % for fields in the
    cols=col0:colL;     % subdomain considered
    
    [li,~]=SPLengthSing(...               % Calculating length of subdomain
        type(i),y(i),gm(i),gm(i+1),Sm(i),Sm(i+1),DX0i,T,tm,Eb);
    
    % Calculating full fields for subdomain i
    [x(:,cols),DX(:,cols),S(:,cols),g(:,cols)]=...
        SPFieldsSing(type(i),Gm(i),y(i),g0i,S0i,DX0i,li,T,tm,Eb,Nx,Ng);
    x(:,cols)=x(:,cols)+repmat(x0i,1,Nx+1);   

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fields for subdomain at the edge of the overlap %%%%%%%%%%%%%%%%%%%%%%%%%

i=sk(n+1);          % Identifying subdomain at the edge

x0i=x(:,colL);      % Calculating fields at
DX0i=DX(:,colL);    % x=x0(i) (which corresponds to xL(i-1))

col0=colL+1;        % Updating column
colL=col0+Nx;       % for fields in the
cols=col0:colL;     % edge subdomain

li=L-x0i;           % Calculating length of subdomain at the edge

if i==N+1           % if A CRACK HAS BEEN FORMED at the edge, use 
                    % expressions for fully debonded subdomain
                    % (Table 1 Eq.h)
                    
    % Calculating full fields for subdomain at the edge
    x(:,cols)=SPlinspace2(0,li,Nx+1);
    DX(:,cols)=repmat(DX0i,1,Nx+1);
    S(:,cols)=0;
    g(:,cols)=gm(end)+DX(:,cols)/(tm*Eb).*x(:,cols);

else                % if NO CRACK has been formed at the edge, 
                    % use general case
    g0i=gm(i)*ones(Ng+1,1);             % Calculating fields at 
    S0i=interp1(gm,Sm,g0i,'linear');    % x=x0(i)
    
    % Calculating full fields for subdomain at the edge
    [x(:,cols),DX(:,cols),S(:,cols),g(:,cols)]=...
    SPFieldsSing(type(i),Gm(i),y(i),g0i,S0i,DX0i,li,T,tm,Eb,Nx,Ng);
end

% Updating x=xi+x0i;
x(:,cols)=x(:,cols)+repmat(x0i,1,Nx+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of fields for edge subdomain            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hence end of fields for multiple subdomains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end