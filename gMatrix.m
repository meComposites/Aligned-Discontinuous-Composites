function [gcritMax] = gMatrix(N,type,Y,gm,Sm,ET,DT,tm,EE,DE)
% Calculates matrix of Length of Process Zones (Eq. 11)
%
% Example for N=4:
%
%          --                               --
%          |  g[11]     g[21]      g[321]    |
% g1[pz] = |  g[12]     g[212]     g[3212]   |
%          |  g[123]    g[2123]    g[32123]  |
%          --                               --
%
% g2 ...
% ...
%
% Vector of active subdomains is s={i, i+1, ... , i+n=j}, where:
% i is the subdomain at the centre of the overlap;
% i+n=j is the subdomain at the edge of the overlap;
% n+1 is the number of active subdomains.
%

[L,~,~,~] = geometry();
opts = optimset('TolX',eps*L/2/10000,'Display','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initiating matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gcritMax=NaN(N,N,N);

for a=1:N
    for b=a:N
        for c=a:N
            
            if a==1
                gini=0;
            else
                gini=gm(a-1);
            end
                
%                 [length_act,jkl,~]= fminsearch(@find_gcrit,L/2,opts);
                [gcrit,jkl,~]= fminbnd(@find_gcrit,gini,gm(a),opts);
                
                %[gcrit,jkl,~]= fzero(@find_gcrit_sgn,[gini gm(a)],opts);
                
                gcritMax(a,b,c)=gcrit;
                
        end
    end
    
    % add shear strain for deactivation (if delamination)
    gcritMax(a,N+1,:)=gm(a);
    gcritMax(a,:,N+1)=gm(a);
    
end


    function find_gcrit = find_gcrit(gcrit_guess)
        
        % Right
        lengthR=0;
        P0=0;
        for j=a:b
            if j==a
                g0=gcrit_guess;
                S0=interp1([0 gm'],[0 Sm'],gcrit_guess,'linear');
            else
                g0=gm(j-1);
                S0=Sm(j-1);
            end
            
            [length,P0]=LengthSing(type(j),Y(j),g0,gm(j),S0,Sm(j),abs(P0),ET,DT,tm,EE,DE);
            
            lengthR=lengthR+length;
        end
        
        % Left
        lengthL=0;
        P0=0;
        for j=a:c
            if j==a
                g0=gcrit_guess;
                S0=interp1([0 gm'],[0 Sm'],gcrit_guess,'linear');
            else
                g0=gm(j-1);
                S0=Sm(j-1);
            end
            
            [length,P0]=LengthSing(type(j),Y(j),g0,gm(j),S0,Sm(j),abs(P0),ET,DT,tm,EE,DE);
            
            lengthL=lengthL+length;
        end
        
        % Error function
        find_gcrit = abs(L/2 - lengthL - lengthR);
        
    end



    function find_gcrit = find_gcrit_sgn(gcrit_guess)
        
        % Right
        lengthR=0;
        P0=0;
        for j=a:b
            if j==1
                g0=gcrit_guess;
                S0=interp1([0 gm'],[0 Sm'],gcrit_guess,'linear');
            else
                g0=gm(j-1);
                S0=Sm(j-1);
            end
            
            [length,P0]=LengthSing(type(j),Y(j),g0,gm(j),S0,Sm(j),abs(P0),ET,DT,tm,EE,DE);
            
            lengthR=lengthR+length;
        end
        
        % Left
        lengthL=0;
        P0=0;
        for j=a:c
            if j==1
                g0=gcrit_guess;
                S0=interp1([0 gm'],[0 Sm'],gcrit_guess,'linear');
            else
                g0=gm(j-1);
                S0=Sm(j-1);
            end
            
            [length,P0]=LengthSing(type(j),Y(j),g0,gm(j),S0,Sm(j),abs(P0),ET,DT,tm,EE,DE);
            
            lengthL=lengthL+length;
        end
        
        % Error function
        find_gcrit = L/2 - lengthL;
        
    end

end