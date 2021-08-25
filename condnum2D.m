clear all
startup; 
curpath = pwd;
addpath(sprintf('%s/%s',curpath,'BF/2D/MBF/src'))
dirs = {'BF/2D/MBF'};
for s = dirs
  addpath(sprintf('%s/%s/src',curpath,s{:}))
  addpath(sprintf('%s/%s/test',curpath,s{:}))
end

addpath(sprintf('%s/%s',curpath,'BF/2D/GBF/src'))
dirs = {'BF/2D/GBF'};
for s = dirs
  addpath(sprintf('%s/%s/src',curpath,s{:}))
  addpath(sprintf('%s/%s/test',curpath,s{:}))
end

funs = {@fun2 @fun0_2d_decay};
func_names = ["fun2", "fun0_2d_decay"];
OutPutFile = fopen(["condnum2D.txt"],'w');


mR = 20;
maxRank = [6 6 6 6];
occ = 4*4;
tol_bf = 1E-6;
tol_peel = 1E-3;
tol_HIF = 1E-3;
maxit1 = 40;
maxit2 = 20;

%dims = 2.^[8 9 10]
dims = [32 64]
cases = length(dims);

% fileID = fopen('meaningless.txt','w');

for ind = 1:length(func_names)
    func_name = func_names(ind);
    fun = funs{ind};
    fprintf(OutPutFile, "\n");
    fprintf(OutPutFile, "\n");
    fprintf(OutPutFile, "kernel: "+func_name+" \n");

    for i = 1:cases
        N = dims(i);
        kN = log2(N);

        fileID = fopen(['results_2d/MBF_',num2str(N^2),'.txt'],'w');

        Hr = 75*(2^kN);
        
        kbox = [-N/2,N/2;-N/2,N/2];
        k = -N/2:N/2-1;
        [k1s,k2s] = ndgrid(k);
        k1s = k1s(:);  k2s = k2s(:);
        kk = [k1s k2s];

        xbox = [0,1;0,1];
        x = (0:N-1)'/N;
        [x1s,x2s] = ndgrid(x);
        x1s = x1s(:);  x2s = x2s(:);
        xx = [x1s x2s];

        fprintf(OutPutFile, 'N: %d \n', N^2);
        Factor = mbf_explicit(fun, xx, xbox, kk, kbox, mR, tol_bf, 1);
        E = eye(N^2);
        A = apply_mbf(Factor, E);
        ATA = apply_mbf_adj(Factor, apply_mbf(Factor, E));
        condA = cond(A);
        condATA = cond(ATA);
        fprintf(OutPutFile, 'condition number of A: %10.4e \n', condA);
        fprintf(OutPutFile, 'condition number of ATA: %10.4e \n', condATA);
    end
end
