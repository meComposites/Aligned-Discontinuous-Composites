function [XeRemoteAll,PropertiesAll,Giicm,XePointsAll,Sgm,... % Always available
    fieldx,fieldg,fieldS,fieldXB,ijActive,gCrit] = SPSL(Ei,Ti,tmi,Lall,SmIn,Giicm,G)   % Only for single L 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[gm,Sm,Giicm,T,tm,Eb,Nx,Ng] = SPInputs(Ei,Ti,tmi,SmIn,Giicm,G);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Preliminary calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,type,Gm,y,gm,Sm,Sgm] = SPPreCalc(gm,Sm,Giicm,T,tm,Eb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Matrix of Length of Process Zones %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lpz = SPLpzMatrix(N,type,y,gm,Sm,T,tm,Eb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Loop for all L required %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(Lall)>1
    Nx=1;       % If there are several L, it will not output fields.
end

[x,DX,S,g,nMax,XeRemoteAll,PropertiesAll,XePointsAll,ijActive,DXL,gCrit] = ...
SPSLLoop(N,type,Gm,y,gm,Sm,Lpz,Lall,T,tm,Eb,Nx,Ng);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Calculates fields if there is only on value of L %%%%%%%%%%%%%%%%%%%%

if length(Lall)==1
    [fieldx,fieldg,fieldS,fieldXB] = ...
        SPPostProcFields(Nx,nMax,x,DX,S,g,DXL);
end

end