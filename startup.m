% add core paths
curpath = pwd;
% dirs = {'Peeling','hifie','rskelf','rskelf/sv','rskelf/mv',...
%     'ButterflyLab-master/v2/matvec/src','IBF/Hmatrix','IBF/Hmatrix_inv',...
%     'HODLR/Construction','HODLR/Inversion'};
dirs = {'Peeling','hifie','rskelf','rskelf/sv','rskelf/mv',...
    'BF/1D/HSS', 'BF/1D/src', 'BF/1D/test', 'BF/2D', 'HQR', 'regularization',...
    'HODLR/Construction','HODLR/Inversion'};
for s = dirs
  addpath(sprintf('%s/%s',curpath,s{:}))
end

clear

