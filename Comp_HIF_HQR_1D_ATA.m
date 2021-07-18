clear all
startup; 
curpath = pwd;
dirs = {'BF/1D'};
for s = dirs
  addpath(sprintf('%s/%s/src',curpath,s{:}));
  addpath(sprintf('%s/%s/test',curpath,s{:}));
end


func_name = '2';%'fun_FIO_var2';'fun_FIO';'fun_FIO_5';'fun_FIO_var4';
OutPutFile = fopen(['comp_1d/Comp_HIFvsHQR_ATA_',func_name,'.txt'],'w');

mR = 24;
occ = 32;
tol_bf = 1E-8;
tol_peel = 1E-8;
tol_RSS = 1E-8;
maxit = 120;
restart1 = 10;
maxit1 = 12;
repeat_num = 1;

% dims = 2.^[8 9 10]
dims = 2.^[8 9 10 11 12 13 14 15 16 17 18 19]
cases = length(dims);
bftime = zeros(cases, 1);
bferr = zeros(cases, 1);
condAs = zeros(cases, 1);
condATAs = zeros(cases, 1);

factime_hodlr = zeros(cases, 1);
facerr_hodlr = zeros(cases, 1);

factime_hif = zeros(cases, 1);
factime_hqr = zeros(cases, 1);
apptime_hif = zeros(cases, 1);
apptime_hqr = zeros(cases, 1);
soltime_hif = zeros(cases, 1);
soltime_hqr = zeros(cases, 1);
apperr_hif = zeros(cases, 1);
apperr_hqr = zeros(cases, 1);
solerr_hif = zeros(cases, 1);
solerr_hqr = zeros(cases, 1);

ranks_hif = zeros(cases, 1);
ranks_hqr = zeros(cases, 1);
solerrpcg = zeros(cases, 1);
solerrpcg_hif = zeros(cases, 1);
solerrpcg_hqr = zeros(cases, 1);
iters = zeros(cases, 1);
iters_hif = zeros(cases, 1);
iters_hqr = zeros(cases, 1);

for i = 1:cases
    N = dims(i)+1;

    fileID = fopen(['results_1d/BFF_ATA_',func_name,'_N_',num2str(N),'_v12.txt'],'w');
    fprintf(fileID,'\n');
    fprintf(fileID,'\n');
    fprintf(fileID,'------------------------------------------\n');
    fprintf(fileID,'N: %10d \n',N-1);
    
    fprintf(OutPutFile, '\n');
    fprintf(OutPutFile, '\n');
    fprintf(OutPutFile, '------------------------------------------\n');
    fprintf(OutPutFile, 'N: %10d \n',N-1);

    rk = mR;
    num = str2double(func_name);
    generate_Z;
    Afunc = @(i,j) Z0fun(i,j)/nZ;
    %% COnstruct BF factorization
    [bftime(i), Factor, bferr(i), FL, FU] = run_bf_explicit_mod(N, Afunc, mR, tol_bf, fileID);
    fprintf(OutPutFile, 'BF Fac err/time: %10.4e/%10.4e (s) \n', bferr(i), bftime(i));
    
    Afun = @(x) LUBF_sol(FL, HSSBF_apply(Factor,LUBF_sol(FU, x, 'U')), 'L');
    ATfun = @(x) LUBF_adj_sol(FU, HSSBF_adj_apply(Factor,LUBF_adj_sol(FL, x, 'L')), 'U');

    %% Construct HODLR factorization using peeling algorithm
    tStart=tic;
    for j = 1:repeat_num
      [F,HODLR] = HODLR_construction( N, Afun, ATfun, tol_peel, fileID, occ, 64,64);
    end
    t = toc(tStart)/repeat_num;
    factime_hodlr(i) = t;
    fprintf(fileID,'time H_matrix construction 10e-6: %10.4e )\n',t);
    lvls = floor(log2(length(HODLR))-1);
    HODLR = HODLR_transfer(HODLR, lvls, lvls, 1); %% use another way to store HODLR
    % [e,niter] = snorm(N,@(x)(HSSBF_adj_apply(Factor,HSSBF_apply(Factor,x)) - hodlr_apply(HODLR,x)),[],[],32);
    % e = e/snorm(N,@(x)(HSSBF_adj_apply(Factor,HSSBF_apply(Factor,x))),[],[],1);
    v = randn(N,1) + 1i*randn(N,1);
    e = norm( ATfun(Afun(v)) - hodlr_apply(HODLR,v))/norm(ATfun(Afun(v)));
    fprintf(fileID,'mv: %10.4e \n',e);
    fprintf(OutPutFile, 'HODLR Fac err/time: %10.4e/%10.4e (s) \n', e, t);
    facerr_hodlr(i) = e;

    f = randn(N,1) + 1i*randn(N,1);
    
    % Construct RSS factorization of HODLR matrix
    tStart_HIF = tic;
    for j = 1:repeat_num
      [G] = RSS_ldl(F,tol_RSS,fileID);
    end
    t = toc(tStart_HIF)/repeat_num;
    factime_hif(i) = t;
    fprintf(fileID,'time HIF construction %10.4e \n',t);
    fprintf(OutPutFile, 'RSS Fac time: %10.4e \n', t);
   
    tic;
    for j = 1:repeat_num
      RSS_apply_lu(G,f);
    end
    t=toc/repeat_num;
    % [e,niter] = snorm(N,@(x)(HSSBF_adj_apply(Factor,HSSBF_apply(Factor,x)) - RSS_apply(G,x)),[],[],32);
    % e = e/snorm(N,@(x)(HSSBF_adj_apply(Factor,HSSBF_apply(Factor,x))),[],[],1);
    e = norm(ATfun(Afun(f)) - RSS_apply_lu(G,f))/norm(ATfun(Afun(f)));
    fprintf(fileID,'mv: %10.4e  time %10.4e\n',e,t);
    apptime_hif(i) = t;
    apperr_hif(i) = e;
    fprintf(OutPutFile, 'RSS mv err/time: %10.4e/%10.4e \n', e, t);

    tic;
    for j = 1:repeat_num
      RSS_inv_lu(G,f);
    end
    t=toc;
    % NORM(INV(A) - INV(F))/NORM(INV(A)) <= NORM(I - A*INV(F))
    %[e,niter] = snorm(N,@(x)(x-RSS_inv(G,HSSBF_adj_apply(Factor,HSSBF_apply(Factor,x)))),[],[],1);
    %fprintf(fileID,'sv: %10.4e / (%4d) time %10.4e\n',e,niter,t);
    %NORM(INV(L*L') - INV(F))/NORM(INV(L*L')) <= NORM(I - L*INV(F)*L')
    % [e,niter] = snorm(N,@(x)(x-HSSBF_apply(Factor,RSS_inv(G,HSSBF_adj_apply(Factor,x)))),[],[],32);
    e = norm(f-Afun(RSS_inv_lu(G,ATfun(f))))/norm(f);
    fprintf(fileID,'sv: %10.4e time %10.4e\n',e,t);
    soltime_hif(i) = t;
    solerr_hif(i) = e;
    fprintf(OutPutFile, 'RSS sv err/time: %10.4e/%10.4e \n', e, t);

    % Construct HODLR-QR factorization of HODLR matrix
    tStart_HQR = tic;
    for j = 1:repeat_num
      [Y, YB, YC, T, R, rk] = hodlrqr(HODLR, [], [], [], lvls, 1, tol_RSS);
    end
    t = toc(tStart_HQR)/repeat_num;
    factime_hqr(i) = t;
    ranks_hqr(i) = rk;
    fprintf(fileID,'time/rank HODLR-QR construction %10.4e/%d \n',t, rk);
    fprintf(OutPutFile, 'HQR Fac time/rank: %10.4e/%d \n', t, rk);
    
    tic;
    for j = 1:repeat_num
      y = hodlrqr_apply(Y, T, R, f);
    end
    t=toc/repeat_num;
    % [e,niter] = snorm(N,@(x)(HSSBF_adj_apply(Factor,HSSBF_apply(Factor,x)) - hodlrqr_apply(Y, T, R, x)),[],[],32);
    % e = e/snorm(N,@(x)(HSSBF_adj_apply(Factor,HSSBF_apply(Factor,x))),[],[],1);
    e = norm(ATfun(Afun(f)) - hodlrqr_apply(Y, T, R, f))/norm(ATfun(Afun(f)));
    fprintf(fileID,'mv: %10.4e time %10.4e\n',e,t);
    apptime_hqr(i) = t;
    apperr_hqr(i) = e;
    fprintf(OutPutFile, 'HQR mv err/time: %10.4e/%10.4e \n', e, t);

    tic;
    for j = 1:repeat_num
      hodlrqr_inv(Y, T, R, f);
    end
    t=toc/repeat_num;
    % [e,niter] = snorm(N,@(x)(x-HSSBF_apply(Factor,hodlrqr_inv(Y, T, R, HSSBF_adj_apply(Factor,x)))),[],[],32);
    e = norm(f-Afun(hodlrqr_inv(Y, T, R, ATfun(f))))/norm(f);
    fprintf(fileID,'sv: %10.4e time %10.4e\n',e,t);
    soltime_hqr(i) = t;
    solerr_hqr(i) = e;
    fprintf(OutPutFile, 'HQR sv err/time: %10.4e/%10.4e \n', e, t);

    %cond(RSS_inv(G, eye(N)))
    %cond(hodlrqr_inv(Y, T, R, eye(N)))
    %Q = eye(N) - hodlr_apply(Y, hodlr_apply(T, hodlr_adj_apply(Y, eye(N))));
    %norm(Q'*Q-eye(N))/norm(eye(N))    

    % run CG 
    for tol=[1E-10,1E-14]
        tic
        [x,flag,relres,iter] = gmres(@(x) Afun(x),LUBF_sol(FL,HSSBF_apply(Factor,f),'L'),restart1,tol,maxit1);
        t1 = toc;
        relerr = norm(f-LUBF_sol(FU, x, 'U'))/norm(f);
        fprintf(fileID,'GMRES without preconditioning: tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t1,iter(1)*iter(2), relerr);
        fprintf(OutPutFile, 'GMRES with no preconditioning: tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t1,iter(1)*iter(2), relerr);
        solerrpcg(i) = relerr;
        iters(i) = iter(1)*iter(2);


        y = HSSBF_apply(Factor,f);
        b = ATfun(LUBF_sol(FL, y, 'L'));

        tic
        [x,flag,relres,iter] = pcg(@(x) ATfun(Afun(x)),b,tol,maxit);
        t1 = toc;
        relerr = norm(f-LUBF_sol(FU, x, 'U'))/norm(f);
        fprintf(fileID,'CG with A* preconditioning: tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t1,iter, relerr);
        fprintf(OutPutFile, 'CG with A*    preconditioning: tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t1,iter, relerr);
        solerrpcg(i) = relerr;
        iters(i) = iter;

        tic
        [x,flag,relres,iter] = pcg(@(x) ATfun(Afun(x)),b,tol,maxit,@(x) RSS_inv_lu(G,x)); 
        t = toc;
        relerr = norm(f-LUBF_sol(FU, x, 'U'))/norm(f);
        fprintf(fileID,'CG with HIF preconditioning: tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t,iter, relerr);
        fprintf(OutPutFile, 'CG with HIF   preconditioning: tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t,iter, relerr);
        solerrpcg_hif(i) = relerr;
        iters_hif(i) = iter;

        tic
        [x,flag,relres,iter] = pcg(@(x) ATfun(Afun(x)),b,tol,maxit,@(x) hodlrqr_inv(Y, T, R, x)); 
        t = toc;
        relerr = norm(f-LUBF_sol(FU, x, 'U'))/norm(f);
        fprintf(fileID,'CG with HQR preconditioning: tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t,iter, relerr);
        fprintf(OutPutFile, 'CG with HQR   preconditioning: tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t,iter, relerr);
        solerrpcg_hqr(i) = relerr;
        iters_hqr(i) = iter;
    end

end
