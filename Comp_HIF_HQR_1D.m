clear all
startup; 
curpath = pwd;
dirs = {'BF/1D'};
for s = dirs
  addpath(sprintf('%s/%s/src',curpath,s{:}));
  addpath(sprintf('%s/%s/test',curpath,s{:}));
end


func_name = 'fun_FIO_var4';%'fun_FIO_var2';'fun_FIO';'fun_FIO_5';'fun_FIO_var4';
OutPutFile = fopen(['comp_1d/Comp_HIFvsHQR_',func_name,'_v12.txt'],'w');

mR = 16;
occ = 32;
tol_bf = 1E-13;
tol_peel = 1E-12;
tol_RSS = 1E-12;
maxit = 120;
repeat_num = 1;

%dims = 2.^[8 9 10]
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
    N = dims(i);

    fileID = fopen(['results_1d/BFF_',func_name,'_N_',num2str(N),'_v12.txt'],'w');
    fprintf(fileID,'\n');
    fprintf(fileID,'\n');
    fprintf(fileID,'------------------------------------------\n');
    fprintf(fileID,'N: %10d \n',N);
    
    fprintf(OutPutFile, '\n');
    fprintf(OutPutFile, '\n');
    fprintf(OutPutFile, '------------------------------------------\n');
    fprintf(OutPutFile, 'N: %10d \n',N);

    %% COnstruct BF factorization
    [bftime(i), Factor, bferr(i)] = run_bf_explicit(N, func_name, mR, tol_bf, fileID);
    fprintf(OutPutFile, 'BF Fac err/time: %10.4e/%10.4e (s) \n', bferr(i), bftime(i));

    %% Construct HODLR factorization using peeling algorithm
    tStart=tic;
    for j = 1:repeat_num
      [F,HODLR] = HODLR_construction( N, @(x) apply_bf(Factor,x), @(x) apply_bf_adj(Factor,x), tol_peel, fileID, occ, 64,64);
    end
    t = toc(tStart)/repeat_num;
    factime_hodlr(i) = t;
    fprintf(fileID,'time H_matrix construction 10e-6: %10.4e )\n',t);
    lvls = floor(log2(length(HODLR))-1);
    HODLR = HODLR_transfer(HODLR, lvls, lvls, 1); %% use another way to store HODLR
    % [e,niter] = snorm(N,@(x)(apply_bf_adj(Factor,apply_bf(Factor,x)) - hodlr_apply(HODLR,x)),[],[],32);
    % e = e/snorm(N,@(x)(apply_bf_adj(Factor,apply_bf(Factor,x))),[],[],1);
    v = randn(N,1) + 1i*randn(N,1);
    e = norm(apply_bf_adj(Factor,apply_bf(Factor,v)) - hodlr_apply(HODLR,v))/norm(apply_bf_adj(Factor,apply_bf(Factor,v)));
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
    % [e,niter] = snorm(N,@(x)(apply_bf_adj(Factor,apply_bf(Factor,x)) - RSS_apply(G,x)),[],[],32);
    % e = e/snorm(N,@(x)(apply_bf_adj(Factor,apply_bf(Factor,x))),[],[],1);
    e = norm(apply_bf_adj(Factor,apply_bf(Factor,f)) - RSS_apply_lu(G,f))/norm(apply_bf_adj(Factor,apply_bf(Factor,f)));
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
    %[e,niter] = snorm(N,@(x)(x-RSS_inv(G,apply_bf_adj(Factor,apply_bf(Factor,x)))),[],[],1);
    %fprintf(fileID,'sv: %10.4e / (%4d) time %10.4e\n',e,niter,t);
    %NORM(INV(L*L') - INV(F))/NORM(INV(L*L')) <= NORM(I - L*INV(F)*L')
    % [e,niter] = snorm(N,@(x)(x-apply_bf(Factor,RSS_inv(G,apply_bf_adj(Factor,x)))),[],[],32);
    e = norm(f-apply_bf(Factor,RSS_inv_lu(G,apply_bf_adj(Factor,f))))/norm(f);
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
    % [e,niter] = snorm(N,@(x)(apply_bf_adj(Factor,apply_bf(Factor,x)) - hodlrqr_apply(Y, T, R, x)),[],[],32);
    % e = e/snorm(N,@(x)(apply_bf_adj(Factor,apply_bf(Factor,x))),[],[],1);
    e = norm(apply_bf_adj(Factor,apply_bf(Factor,f)) - hodlrqr_apply(Y, T, R, f))/norm(apply_bf_adj(Factor,apply_bf(Factor,f)));
    fprintf(fileID,'mv: %10.4e time %10.4e\n',e,t);
    apptime_hqr(i) = t;
    apperr_hqr(i) = e;
    fprintf(OutPutFile, 'HQR mv err/time: %10.4e/%10.4e \n', e, t);

    tic;
    for j = 1:repeat_num
      hodlrqr_inv(Y, T, R, f);
    end
    t=toc/repeat_num;
    % [e,niter] = snorm(N,@(x)(x-apply_bf(Factor,hodlrqr_inv(Y, T, R, apply_bf_adj(Factor,x)))),[],[],32);
    e = norm(f-apply_bf(Factor,hodlrqr_inv(Y, T, R, apply_bf_adj(Factor,f))))/norm(f);
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
        b = apply_bf_adj(Factor,apply_bf(Factor,f));

        tic
        [x,flag,relres,iter] = pcg(@(x) apply_bf_adj(Factor,apply_bf(Factor,x)),b,tol,maxit);
        t1 = toc;
        relerr = norm(f-x)/norm(f);
        fprintf(fileID,'CG with A* preconditioning: tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t1,iter, relerr);
        fprintf(OutPutFile, 'CG with A* preconditioning : tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t1,iter, relerr);
        solerrpcg(i) = relerr;
        iters(i) = iter;

        tic
        [x,flag,relres,iter] = pcg(@(x) apply_bf_adj(Factor,apply_bf(Factor,x)),b,tol,maxit,@(x) RSS_inv_lu(G,x)); 
        t = toc;
        relerr = norm(f-x)/norm(f);
        fprintf(fileID,'CG with HIF preconditioning: tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t,iter, relerr);
        fprintf(OutPutFile, 'CG with HIF preconditioning: tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t,iter, relerr);
        solerrpcg_hif(i) = relerr;
        iters_hif(i) = iter;

        tic
        [x,flag,relres,iter] = pcg(@(x) apply_bf_adj(Factor,apply_bf(Factor,x)),b,tol,maxit,@(x) hodlrqr_inv(Y, T, R, x)); 
        t = toc;
        relerr = norm(f-x)/norm(f);
        fprintf(fileID,'CG with HQR preconditioning: tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t,iter, relerr);
        fprintf(OutPutFile, 'CG with HQR preconditioning: tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t,iter, relerr);
        solerrpcg_hqr(i) = relerr;
        iters_hqr(i) = iter;
    end

end


N = dims(2:end);
logN = log2(N);
NlogN = logN + log2(logN);
N2logN = logN + 2*log2(logN);
N4logN = logN + 4*log2(logN);

bftime = bftime(2:end);
bferr = bferr(2:end);
factime_hodlr = factime_hodlr(2:end);
facerr_hodlr = facerr_hodlr(2:end);
factime_hif = factime_hif(2:end);
factime_hqr = factime_hqr(2:end);
apptime_hif = apptime_hif(2:end);
apptime_hqr = apptime_hqr(2:end);
soltime_hif = soltime_hif(2:end);
soltime_hqr = soltime_hqr(2:end);
apperr_hif = apperr_hif(2:end);
apperr_hqr = apperr_hqr(2:end);
solerr_hif = solerr_hif(2:end);
solerr_hqr = solerr_hqr(2:end);

ranks_hif = ranks_hif(2:end);
ranks_hqr = ranks_hqr(2:end);
solerrpcg = solerrpcg(2:end);
solerrpcg_hif = solerrpcg_hif(2:end);
solerrpcg_hqr = solerrpcg_hqr(2:end);


fig = figure(1);
hold on;
h(1) = plot(logN, log2(bftime));
h(2) = plot(logN, NlogN-NlogN(1)+log2(bftime(1)));
h(3) = plot(logN, N2logN-N2logN(1)+log2(bftime(1)));
xlabel('Log(N)');
ylabel('Log(Time)/s'); 
title('BF time scaling');
legend(h, 'BF time', 'N log N', 'N log^2 N');
% legend({'BF time', 'N log N', 'N log^2 N'}, 'Location', 'northwest');
axis square;
saveas(fig, "./comp/1D/bftime_" + func_name + ".png");
hold off;

fig = figure(2);
hold on;
h(1) = plot(logN, log10(bferr));
xlabel('Log(N)');
ylabel('Log10(error)'); 
title('BF error');
% legend({'BF time', 'N log N', 'N log^2 N'}, 'Location', 'northwest');
axis square;
saveas(fig, "./comp/1D/bferr_" + func_name + ".png");
hold off;


fig = figure(3);
hold on;
h(1) = plot(logN, log2(factime_hodlr));
h(2) = plot(logN, NlogN-NlogN(1)+log2(factime_hodlr(1)));
h(3) = plot(logN, N2logN-N2logN(1)+log2(factime_hodlr(1)));
h(4) = plot(logN, N4logN-N4logN(1)+log2(factime_hodlr(1)));
xlabel('Log(N)');
ylabel('Log(Time)/s'); 
title('Peeling (HODLR) time scaling');
legend({'Peeling', 'N log N', 'N log^2 N', 'N log^4 N'}, 'Location', 'northwest');
axis square;
saveas(fig, "./comp/1D/hodlrtime_" + func_name + ".png");
hold off;

fig = figure(4);
hold on;
h(1) = plot(logN, log10(facerr_hodlr));
xlabel('Log(N)');
ylabel('Log10(Error)'); 
title('Peeling (HODLR) error');
% legend({'Peeling', 'N log N', 'N log^2 N'}, 'Location', 'northwest');
axis square;
saveas(fig, "./comp/1D/hodlrerr_" + func_name + ".png");
hold off;

fig = figure(5);
hold on;
h(1) = plot(logN, log2(factime_hif));
h(2) = plot(logN, log2(factime_hqr));
h(3) = plot(logN, NlogN-NlogN(1)+(log2(factime_hif(1))+log2(factime_hqr(1)))/2);
h(4) = plot(logN, N2logN-N2logN(1)+(log2(factime_hif(1))+log2(factime_hqr(1)))/2);
xlabel('Log(N)');
ylabel('Log(Time)/s'); 
title('Fac time scaling');
legend({'HIF', 'HQR', 'N log N', 'N log^2 N'}, 'Location', 'northwest');
axis square;
saveas(fig, "./comp/1D/factime_" + func_name + ".png");
hold off;

fig = figure(6);
hold on;
h(1) = plot(logN, log2(apptime_hif));
h(2) = plot(logN, log2(apptime_hqr));
h(3) = plot(logN, NlogN-NlogN(1)+(log2(apptime_hif(1))+log2(apptime_hqr(1)))/2);
h(4) = plot(logN, N2logN-N2logN(1)+(log2(apptime_hif(1))+log2(apptime_hqr(1)))/2);
xlabel('Log(N)');
ylabel('Log(Time)/s'); 
title('App time scaling');
legend({'HIF', 'HQR', 'N log N', 'N log^2 N'}, 'Location', 'northwest');
axis square;
saveas(fig, "./comp/1D/apptime_" + func_name + ".png");
hold off;

fig = figure(7);
hold on;
h(1) = plot(logN, log2(soltime_hif));
h(2) = plot(logN, log2(soltime_hqr));
h(3) = plot(logN, NlogN-NlogN(1)+(log2(soltime_hif(1))+log2(soltime_hqr(1)))/2);
h(4) = plot(logN, N2logN-N2logN(1)+(log2(soltime_hif(1))+log2(soltime_hqr(1)))/2);
xlabel('Log(N)');
ylabel('Log(Time)/s'); 
title('Sol time scaling');
legend({'HIF', 'HQR', 'N log N', 'N log^2 N'}, 'Location', 'northwest');
axis square;
saveas(fig, "./comp/1D/soltime_" + func_name + ".png");
hold off;

fig = figure(8);
hold on;
h(1) = plot(logN, log10(apperr_hif));
h(2) = plot(logN, log10(apperr_hqr));
xlabel('Log(N)');
ylabel('Log10(Error)'); 
title('App error');
legend({'HIF', 'HQR'}, 'Location', 'northwest');
axis square;
saveas(fig, "./comp/1D/apperr_" + func_name + ".png");
hold off;

fig = figure(9);
hold on;
h(1) = plot(logN, log10(solerr_hif));
h(2) = plot(logN, log10(solerr_hqr));
xlabel('Log(N)');
ylabel('Log10(Error)'); 
title('Sol error');
legend({'HIF', 'HQR'}, 'Location', 'northwest');
axis square;
saveas(fig, "./comp/1D/solerr_" + func_name + ".png");
hold off;

exit;


