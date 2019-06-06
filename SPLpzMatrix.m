function Lpz = SPLpzMatrix(N,type,y,gm,Sm,T,tm,Eb)
% Calculates matrix of Length of Process Zones (Eq. 11)
%
% Example for N=4:
%
%          --                                    --
%          |  lpz(2)     lpz(2+3)     lpz(2+3+4)  |
% L[pz] =  |  lpz(3)     lpz(3+4)                 |
%          |  lpz(4)                              |
%          --                                    --
%
% Vector of active subdomains is s={i, i+1, ... , i+n=j}, where:
% i is the subdomain at the centre of the overlap;
% i+n=j is the subdomain at the edge of the overlap;
% n+1 is the number of active subdomains.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initiating matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lpz=NaN(N-1,N-1);
DX0pz=zeros(N-1,N-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations for central region only, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n=0 (1st column)                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for line=1:N-1               
    i=line+1;           
    [Lpz(line,1),DX0pz(line,1)]=...
        SPLengthSingCentre(type(i),y(i),gm(i),gm(i+1),Sm(i),Sm(i+1),T,tm,Eb);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-central regions, n>=1                                   %%%%%%%%%%%%%
% (n_max=N-1, when N-1 subdomains (excluding i=1) are active) %%%%%%%%%%%%%
for col=2:N-1
    n=col-1;
    for line=1:(N-1)-n      % total of N-1 subdomains to be grouped in 
                            % sets of n+1 active subdomains
        j=line+col;
        [Lpz(line,col),DX0pz(line,col)]=...
            SPLengthSing(type(j),y(j),gm(j),gm(j+1),Sm(j),Sm(j+1),...
            DX0pz(line,col-1),T,tm,Eb);
        Lpz(line,col)=Lpz(line,col)+Lpz(line,col-1);
    end
end

end