function [S,g,G] = findLastg(MLawS,MLawg,MLawG,Gic,tm)

% MLawSDiff=[MLawS(1) MLawS(2:end)-MLawS(1:end-1)];
% dg=MLawSDiff./MLawG;
% %g=MLawS./MLawG;
% g=cumsum(dg);
g=[0 MLawg];
S=[0 MLawS];

sum=0;

for i=1:length(g)-1
    sum=sum+(S(i+1)+S(i))*(g(i+1)-g(i)); % Double mais divise par deux en dessous
end

g=[g(2:end) g(end)+(2*Gic/tm-sum)/S(end) 10^3]';
S=[S(2:end) 0 0]';
G=[MLawG (S(end-1)-S(end-2))/(g(end-1)-g(end-2)) 0]';

% figure(17)
% hold on
% plot(g,S)

end