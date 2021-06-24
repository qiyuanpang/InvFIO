% add core paths
curpath = pwd;
% dirs = {'Peeling','hifie','rskelf','rskelf/sv','rskelf/mv',...
%     'ButterflyLab-master/v2/matvec/src','IBF/Hmatrix','IBF/Hmatrix_inv',...
%     'HODLR/Construction','HODLR/Inversion'};
dirs = {'Peeling','hifie','rskelf','rskelf/sv','rskelf/mv',...
    'BF', 'HQR',...
    'HODLR/Construction','HODLR/Inversion'};
for s = dirs
  addpath(sprintf('%s/%s',curpath,s{:}))
end

clear

