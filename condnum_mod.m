clear all
startup; 
curpath = pwd;
dirs = {'BF/1D'};
for s = dirs
  addpath(sprintf('%s/%s/src',curpath,s{:}));
  addpath(sprintf('%s/%s/test',curpath,s{:}));
end


func_names = ["fun_FIO_var4","fun_FIO_5","fun_FIO_var1","fun_FIO_var2","fun_FIO"];
OutPutFile = fopen(["condnum_mod.txt"],'w');

mR = 8;
occ = 32;
tol_bf = 1E-12;
tol_peel = 1E-5;
tol_RSS = 1E-4;
maxit = 200;
repeat_num = 1;

%dims = 2.^[8 9 10]
dims = 2.^[8 9 10]
cases = length(dims);

fileID = fopen('meaningless.txt','w');

for func_name = func_names
    fprintf(OutPutFile, "\n");
    fprintf(OutPutFile, "\n");
    fprintf(OutPutFile, "kernel: "+func_name+" \n");

    for i = 1:cases
        N = dims(i);

        fprintf(OutPutFile, 'N: %d \n', N);
        [bftime, Factor, bferr, FL, FU] = run_bf_explicit(N, func_name, mR, tol_bf, fileID);
        E = eye(N);
        % norm(HSSBF_apply(Factor, E) - HSSBF_apply(FL,HSSBF_apply(FU,E)))/norm(HSSBF_apply(Factor, E))
        Afun = @(x) LUBF_sol(FL, HSSBF_apply(Factor,LUBF_sol(FU, x, 'U')), 'L');
        ATfun = @(x) LUBF_adj_sol(FU, HSSBF_adj_apply(Factor,LUBF_adj_sol(FL, x, 'L')), 'U');
        A = Afun(E);
        ATA = ATfun(Afun(E));
        condA = cond(A);
        condATA = cond(ATA);
        fprintf(OutPutFile, 'condition number of A: %10.4e \n', condA);
        fprintf(OutPutFile, 'condition number of ATA: %10.4e \n', condATA);
    end
end