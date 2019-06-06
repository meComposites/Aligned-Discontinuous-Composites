function [g0Crit] = SPTransMult(type,y,gm,Sm,sk,g0prev,L,T,tm,Eb)
% Calculates the transition strain at the centre of the overlap 
% (multiple subdomains in the overlapping length, with transition 
% corresponding to the activation of the next subdomain).
% Based on Table 4, Eq.e.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=length(sk)-1;

g0est=[g0prev,gm(sk(1)+1)]; % Estimated interval for the transition strain:
                            % g0Crit>g0(k) (from previous transition);
                            % g0Crit<gm(subdomain in centre) (corresponds 
                            % to deactivation of the central subdomain);

l=zeros(n+1,1);             % Initiating vector for length of subdomains
DX0=zeros(n+1,1);           % and DXi(xi=0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving for g0Crit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts = optimset('TolX',2*eps(min(g0est)));  % Seting criterion for finding 
                                            % g0Crit numerically

g0Crit = fzero(@Lerror,g0est,opts);         % Finding g0Crit numerically,
                                            % by zeroing
                                            % Lreal-Lcalculated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining function for calculating the error between %%%%%%%%%%%%%%%%%%%%%
% Lreal and Lcalculated                               %%%%%%%%%%%%%%%%%%%%%

    function Lerror = Lerror(g0)
        S0=interp1(gm,Sm,g0,'linear');          % Shear stress at x=0

        [l(1),DX0(1)] = SPLengthSingCentre(...    % Length (partial) of 
            type(sk(1)),y(sk(1)),...            % central subdomain
            g0,gm(sk(1)+1),S0,Sm(sk(1)+1),...
            T,tm,Eb);

        for i=2:n+1
            [l(i),DX0(i)] = SPLengthSing(...      % Full length of 
                type(sk(i)),y(sk(i)),...        % non-central subdomains
                gm(sk(i)),gm(sk(i)+1),...       % (need to find when the
                Sm(sk(i)),Sm(sk(i)+1),...       % last subdomain j becomes
                DX0(i-1),T,tm,Eb);              % complete, from gm(j) 
                                                % to gm(j+1)
        end
        
    Lerror = (sum(l)-L);        % Calculating error of length occupied 
                                % by the n+1 subdomains 
                                % (from g0=g0Crit to g(j+1) )
    end

end