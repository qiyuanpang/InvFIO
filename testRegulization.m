clear all
startup; 
curpath = pwd;
dirs = {'BF/1D'};
for s = dirs
  addpath(sprintf('%s/%s/src',curpath,s{:}));
  addpath(sprintf('%s/%s/test',curpath,s{:}));
end


func_name = 'fun_FIO_var4';%'fun_FIO_var2';'fun_FIO';'fun_FIO_5';'fun_FIO_var4';
OutPutFile = fopen(['comp_1d/Regularization_',func_name,'_L1.txt'],'w');

mR = 20;
occ = 32;
tol_bf = 1E-8;
tol_peel = 1E-5;
tol_RSS = 1E-5;
maxit1 = 40;
maxit2 = 20;
repeat_num = 1;
delta = 10;
mu = 2;
lambda = 2;
regm = 'L1';


%dims = 2.^[8 9 10]
dims = 2.^[8 9 10 11 12 13 14 15 16 17 18]
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

    fileID = fopen(['results_1d/BFF_',func_name,'_N_',num2str(N),'.txt'],'w');
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
      if strcmp(regm, 'L1')
        [F,HODLR] = HODLR_construction( N, @(x)apply_bf(Factor, x), @(x) apply_bf_adj(Factor,x), tol_peel, fileID, occ, 64,64);
      elseif strcmp(regm, 'TV-L1')
        [F,HODLR] = HODLR_construction( N, @(x)apply_bf(Factor, fft(x)), @(x) N*ifft(apply_bf_adj(Factor,x)), tol_peel, fileID, occ, 64,64);
      end
    end
    t = toc(tStart)/repeat_num;
    factime_hodlr(i) = t;
    fprintf(fileID,'time H_matrix construction 10e-6: %10.4e )\n',t);
    lvls = floor(log2(length(HODLR))-1);
    HODLR = HODLR_transfer(HODLR, lvls, lvls, 1); %% use another way to store HODLR
    % [e,niter] = snorm(N,@(x)(apply_bf_adj(Factor,apply_bf(Factor,x)) - hodlr_apply(HODLR,x)),[],[],32);
    % e = e/snorm(N,@(x)(apply_bf_adj(Factor,apply_bf(Factor,x))),[],[],1);
    v = randn(N,1) + 1i*randn(N,1);
    if strcmp(regm, 'L1')
      e = norm(apply_bf_adj(Factor, apply_bf(Factor,v)) - hodlr_apply(HODLR,v))/norm(apply_bf_adj(Factor, apply_bf(Factor,v)));
    elseif strcmp(regm, 'TV-L1')
      e = norm(N*ifft(apply_bf_adj(Factor, apply_bf(Factor,fft(v)))) - hodlr_apply(HODLR,v))/norm(N*ifft(apply_bf_adj(Factor, apply_bf(Factor,fft(v)))));
    end
    fprintf(fileID,'mv: %10.4e \n',e);
    fprintf(OutPutFile, 'HODLR Fac err/time: %10.4e/%10.4e (s) \n', e, t);
    facerr_hodlr(i) = e;
    
    if strcmp(regm, 'L1')
      sp = zeros(N,1);
      for i = ceil(rand(1,1)*N)
        sp(i) = 1;
      end
      f = randn(N,1) + 1i*randn(N,1);
      f = f.*sp;
    elseif strcmp(regm, 'TV-L1')
      f = ones(N,1);
    end
    % sum(abs(f) ~= 0)/N
    
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
    if strcmp(regm, 'L1')
      e = norm(apply_bf_adj(Factor,apply_bf(Factor,f)) - RSS_apply_lu(G,f))/norm(apply_bf_adj(Factor,apply_bf(Factor,f)));
    elseif strcmp(regm, 'TV-L1')
      e = norm(N*ifft(apply_bf_adj(Factor,apply_bf(Factor,fft(f)))) - RSS_apply_lu(G,f))/norm(N*ifft(apply_bf_adj(Factor,apply_bf(Factor,fft(f)))));
    end
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
    if strcmp(regm, 'L1')
      e = norm(f-apply_bf(Factor,RSS_inv_lu(G,apply_bf_adj(Factor,f))))/norm(f);
    elseif strcmp(regm, 'TV-L1')
      e = norm(f-apply_bf(Factor,fft(RSS_inv_lu(G,N*ifft(apply_bf_adj(Factor,f))))))/norm(f);
    end
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
    if strcmp(regm, 'L1')
      e = norm(apply_bf_adj(Factor,apply_bf(Factor,f)) - hodlrqr_apply(Y, T, R, f))/norm(apply_bf_adj(Factor,apply_bf(Factor,f)));
    elseif strcmp(regm, 'TV-L1')
      e = norm(N*ifft(apply_bf_adj(Factor,apply_bf(Factor,fft(f)))) - hodlrqr_apply(Y, T, R, f))/norm(N*ifft(apply_bf_adj(Factor,apply_bf(Factor,fft(f)))));
    end
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
    if strcmp(regm, 'L1')
      e = norm(f-apply_bf(Factor,hodlrqr_inv(Y, T, R, apply_bf_adj(Factor,f))))/norm(f);
    elseif strcmp(regm, 'TV-L1')
      e = norm(f-apply_bf(Factor,fft(hodlrqr_inv(Y, T, R, N*ifft(apply_bf_adj(Factor,f))))))/norm(f);
    end
    fprintf(fileID,'sv: %10.4e time %10.4e\n',e,t);
    soltime_hqr(i) = t;
    solerr_hqr(i) = e;
    fprintf(OutPutFile, 'HQR sv err/time: %10.4e/%10.4e \n', e, t);

    %cond(RSS_inv(G, eye(N)))
    %cond(hodlrqr_inv(Y, T, R, eye(N)))
    %Q = eye(N) - hodlr_apply(Y, hodlr_apply(T, hodlr_adj_apply(Y, eye(N))));
    %norm(Q'*Q-eye(N))/norm(eye(N))    

    % run CG 
    if strcmp(regm,'TV-L1')
      Afun = @(x) mu*N*ifft(apply_bf_adj(Factor,apply_bf(Factor,fft(x)))) - lambda*laplacianp(x,1/N);
      fprintf(fileID, 'Information for TV-L1 \n')
      [F_A,HODLR_A] = HODLR_construction_hmt( N, Afun, tol_peel, fileID, occ, 64,64);
      [G_A] = RSS_ldl(F_A,tol_RSS,fileID);
      HODLR_A = HODLR_transfer(HODLR_A, lvls, lvls, 1);
      fprintf('TV-L1: %10.4e / %10.4e \n', norm(Afun(f)-RSS_apply_lu(G_A,f))/norm(Afun(f)), norm(Afun(f)-hodlr_apply(HODLR_A,f))/norm(Afun(f)))
      [Y_A, YB_A, YC_A, T_A, R_A, rk_A] = hodlrqr(HODLR_A, [], [], [], lvls, 1, tol_RSS);
      b = apply_bf(Factor,fft(f));
    elseif  strcmp(regm,'L1')
      Afun = @(x) mu*apply_bf_adj(Factor, apply_bf(Factor, x)) + lambda*x;
      fprintf(fileID, 'Information for L1 \n')
      [F_A,HODLR_A] = HODLR_construction_hmt( N, Afun, tol_peel, fileID, occ, 64,64);
      [G_A] = RSS_ldl(F_A,tol_RSS,fileID);
      HODLR_A = HODLR_transfer(HODLR_A, lvls, lvls, 1);
      fprintf('L1: %10.4e / %10.4e \n', norm(Afun(f)-RSS_apply_lu(G_A,f))/norm(Afun(f)), norm(Afun(f)-hodlr_apply(HODLR_A,f))/norm(Afun(f)))
      [Y_A, YB_A, YC_A, T_A, R_A, rk_A] = hodlrqr(HODLR_A, [], [], [], lvls, 1, tol_RSS);
      b = apply_bf(Factor, f);
    end
    for tol=[1E-10,1E-12]

        % tic
        % % [x,flag,relres,iter] = pcg(@(x) apply_bf_adj(Factor,apply_bf(Factor,x)),b,tol,maxit);
        % % [x, iter] = LinearizedBregman(@(x)apply_bf_adj(Factor, x), @(x)apply_bf(Factor, x), b, delta, mu, tol, maxit);
        % [x, iter] = SplitBregman(@(x)apply_bf_adj(Factor, x), @(x)apply_bf(Factor, x), b, @(x)x, @(x)eye(N), lambda, tol, maxit);
        % t1 = toc;
        % relerr = norm(f-x)/norm(f);
        % fprintf(fileID,'Linearized Bregman without intial guess     : tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t1,iter, relerr);
        % fprintf(OutPutFile, 'Linearized Bregman without intial guess: tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t1,iter, relerr);
        % solerrpcg(i) = relerr;
        % iters(i) = iter;

        tic
        % [x,flag,relres,iter] = pcg(@(x) apply_bf_adj(Factor,apply_bf(Factor,x)),b,tol,maxit,@(x) RSS_inv_lu(G,x)); 
        % [x, iter] = LinearizedBregman(@(x)apply_bf_adj(Factor, x), @(x)apply_bf(Factor, x), b, delta, mu, tol, maxit, apply_bf(Factor, RSS_inv_lu(G, b)));
        % norm(apply_bf_adj(Factor, apply_bf(Factor, RSS_inv_lu(G, b)))-b)^2/norm(b)^2
        if strcmp(regm,'TV-L1')
          [x, iter] = SplitBregman(regm, mu, lambda, 1/N, @(x)N*ifft(apply_bf_adj(Factor,x)), @(x)hodlr_apply(HODLR_A, x), b, @(x)gradientx(x, 1/N), maxit1, maxit2, tol, @(x)RSS_inv_lu(G_A, x), RSS_inv_lu(G, N*ifft(apply_bf_adj(Factor, b))), RSS_inv_lu(G, N*ifft(apply_bf_adj(Factor, b))), zeros(N,1));
        elseif strcmp(regm,'L1')
          [x, iter] = SplitBregman(regm, mu, lambda, 1/N, @(x)apply_bf_adj(Factor,x), @(x)hodlr_apply(HODLR_A, x), b, @(x)x, maxit1, maxit2, tol, @(x)RSS_inv_lu(G_A,x), RSS_inv_lu(G, apply_bf_adj(Factor, b)), RSS_inv_lu(G, apply_bf_adj(Factor, b)), zeros(N,1));
        end
        t = toc;
        relerr = norm(f-x)/norm(f);
        fprintf(fileID,'Split Bregman with HIF guess           : tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t,iter, relerr);
        fprintf(OutPutFile, 'Split Bregman with HIF guess      : tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t,iter, relerr);
        solerrpcg_hif(i) = relerr;
        iters_hif(i) = iter;

        tic
        % [x,flag,relres,iter] = pcg(@(x) apply_bf_adj(Factor,apply_bf(Factor,x)),b,tol,maxit,@(x) hodlrqr_inv(Y, T, R, x)); 
        % [x, iter] = LinearizedBregman(@(x)apply_bf_adj(Factor, x), @(x)apply_bf(Factor, x), b, delta, mu, tol, maxit, apply_bf(Factor, hodlrqr_inv(Y, T, R, b)));
        % [x, iter] = SplitBregman(@(x)apply_bf_adj(Factor, x), @(x)apply_bf(Factor, x), b, @(x)x, @(x)eye(N), lambda, tol, maxit, apply_bf(Factor, hodlrqr_inv(Y, T, R, b)), apply_bf(Factor, hodlrqr_inv(Y, T, R, b)));
        if strcmp(regm,'TV-L1')
          [x, iter] = SplitBregman(regm, mu, lambda, 1/N, @(x)N*ifft(apply_bf_adj(Factor,x)), @(x)hodlr_apply(HODLR_A, x), b, @(x)gradientx(x, 1/N), maxit1, maxit2, tol, @(x)hodlrqr_inv(Y_A, T_A, R_A, x), hodlrqr_inv(Y, T, R, N*ifft(apply_bf_adj(Factor, b))), hodlrqr_inv(Y, T, R, N*ifft(apply_bf_adj(Factor, b))), zeros(N,1));
        elseif strcmp(regm,'L1')
          [x, iter] = SplitBregman(regm, mu, lambda, 1/N, @(x)apply_bf_adj(Factor,x), @(x)hodlr_apply(HODLR_A, x), b, @(x)x, maxit1, maxit2, tol, @(x)hodlrqr_inv(Y_A, T_A, R_A, x), hodlrqr_inv(Y, T, R, apply_bf_adj(Factor, b)), hodlrqr_inv(Y, T, R, apply_bf_adj(Factor, b)), zeros(N,1));
        end
        t = toc;
        relerr = norm(f-x)/norm(f);
        fprintf(fileID,'Split Bregman with HQR guess           : tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t,iter, relerr);
        fprintf(OutPutFile, 'Split Bregman with HQR guess      : tol %10.2e,   time/#iter %10.4e / %4d, relerr %10.4e \n',tol,t,iter, relerr);
        solerrpcg_hqr(i) = relerr;
        iters_hqr(i) = iter;
    end

end
