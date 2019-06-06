function [Lpz,P0pz] = LpzMatrix(N,type,Y,gm,Sm,ET,DT,tm,EE,DE)
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
P0pz=NaN(N-1,N-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations for central region only, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n=0 (1st column)                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for line=1:N-1               
    i=line+1;           
    [Lpz(line,1),P0pz(line,1)]=...
        LengthSingCentre(type(i),Y(i),gm(i-1),gm(i),Sm(i-1),Sm(i),ET,DT,tm,EE,DE);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-central regions, n>=1                                   %%%%%%%%%%%%%
% (n_max=N-1, when N-1 subdomains (excluding i=1) are active) %%%%%%%%%%%%%
for col=2:N-1
    n=col-1;
    for line=1:(N-1)-n      % total of N-1 subdomains to be grouped in 
                            % sets of n+1 active subdomains
        j=line+col;
        [Lpz(line,col),P0pz(line,col)]=...
            LengthSing(type(j),Y(j),gm(j-1),gm(j),Sm(j-1),Sm(j),...
            P0pz(line,col-1),ET,DT,tm,EE,DE);
        Lpz(line,col)=Lpz(line,col)+Lpz(line,col-1);
    end
end

end