clear all
startup; 
curpath = pwd;
dirs = {'BF/1D'};
for s = dirs
  addpath(sprintf('%s/%s/src',curpath,s{:}));
  addpath(sprintf('%s/%s/test',curpath,s{:}));
end

initf = 'randn'
regm = 'L1'
whichlambda = 'lmus'

func_name = 'fun_FIO_5';%'fun_FIO_var2';'fun_FIO';'fun_FIO_5';'fun_FIO_var4';
OutPutFile = fopen(['comp_1d/Regularization5_save1_',func_name,'_',regm,'_',initf,'_',whichlambda,'.txt'],'w');


mR = 20;
occ = 32;
tol_bf = 1E-13;
tol_peel = 1E-12;
tol_RSS = 1E-12;
maxit1 = 40;
maxit2 = 20;
repeat_num = 1;
delta = 10;
mus = [2^6];
%mus = [2^(-10) 2^(-9) 2^(-8) 2^(-7) 2^(-6) 2^(-5) 2^(-4) 2^(-3) 2^(-2) 2^(-1) 2^0 2^1 2^2 2^3 2^4 2^5 2^6 2^7 2^8 2^9 2^10];
%mus = linspace(2^4, 2^6, 10);
lambdas = [2^(-16) 2^(-12) 2^(-8) 2^(-4) 2^0 2^4 2^(8) 2^(12) 2^(16)];



%dims = 2.^[8 9 10]
dims = 2.^[14 15 16 17]
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

    
    if strcmp(regm, 'L1')
      f = zeros(N,1);
      f(N/4) = 1;
    elseif strcmp(regm, 'TV-L1')
      f = ones(N,1);
      f(N/2+1:N) = 2;
    elseif strcmp(regm, 'F-L1')
      f1 = zeros(N,1);
      f1(N/4) = 1;
      f = ifft(f1);
    end
    % sum(abs(f) ~= 0)/N
       

    % run SB
    fprintf(OutPutFile, 'regularizer: %s  init-scheme: %s \n', regm, initf)
    for mu = mus
      %lmus = [mu/512 mu/256 mu/128 mu/64 mu/32 mu/16 mu/8 mu/4 mu/2 mu mu*2 mu*4 mu*8 mu*16 mu*32 mu*64];
      lmus = [mu/128];
      if strcmp(whichlambda,'ind')
        LAMBDAS = lambdas;
      elseif strcmp(whichlambda, 'lmus');
        LAMBDAS = lmus;
      end
      for lambda = LAMBDAS
        fprintf(OutPutFile, 'mu: 2^{%4d}  lambda: 2^{%4d} \n', log2(mu), log2(lambda))
        if strcmp(regm,'TV-L1')
          Afun = @(x) (mu*N*ifft(apply_bf_adj(Factor,apply_bf(Factor,fft(x)))) + lambda*laplacianp_rc(x,1/N))/N^2;
          fprintf(fileID, 'Information for TV-L1 \n')
          [F_A,HODLR_A] = HODLR_construction_hmt( N, Afun, tol_peel, fileID, occ, 2000,2000);
          [G_A] = RSS_ldl(F_A,tol_RSS,fileID);
          lvls = floor(log2(length(HODLR_A))-1);
          HODLR_A = HODLR_transfer_hmt(HODLR_A, lvls, lvls, 1);
          fprintf(OutPutFile,'TV-L1: %10.4e / %10.4e \n', norm(Afun(f)-RSS_apply_lu(G_A,f))/norm(Afun(f)), norm(Afun(f)-hodlr_apply(HODLR_A,f))/norm(Afun(f)))
          [Y_A, YB_A, YC_A, T_A, R_A, rk_A] = hodlrqr(HODLR_A, [], [], [], lvls, 1, tol_RSS);
          b = apply_bf(Factor,fft(f))/N;
        elseif strcmp(regm,'F-L1')
          Afun = @(x) mu*N*ifft(apply_bf_adj(Factor,apply_bf(Factor,fft(x)))) + lambda*N*x;
          fprintf(fileID, 'Information for T-L1 \n')
          [F_A,HODLR_A] = HODLR_construction_hmt( N, Afun, tol_peel, fileID, occ, 2000,2000);
          [G_A] = RSS_ldl(F_A,tol_RSS,fileID);
          lvls = floor(log2(length(HODLR_A))-1);
          HODLR_A = HODLR_transfer_hmt(HODLR_A, lvls, lvls, 1);
          fprintf(OutPutFile,'F-L1: %10.4e / %10.4e \n', norm(Afun(f)-RSS_apply_lu(G_A,f))/norm(Afun(f)), norm(Afun(f)-hodlr_apply(HODLR_A,f))/norm(Afun(f)))
          [Y_A, YB_A, YC_A, T_A, R_A, rk_A] = hodlrqr(HODLR_A, [], [], [], lvls, 1, tol_RSS);
          b = apply_bf(Factor,fft(f));
        elseif  strcmp(regm,'L1')
          Afun = @(x) mu*apply_bf_adj(Factor, apply_bf(Factor, x)) + lambda*x;
          fprintf(fileID, 'Information for L1 \n')
          [F_A,HODLR_A] = HODLR_construction_hmt( N, Afun, tol_peel, fileID, occ, 2000,2000);
          [G_A] = RSS_ldl(F_A,tol_RSS,fileID);
          lvls = floor(log2(length(HODLR_A))-1);
          HODLR_A = HODLR_transfer_hmt(HODLR_A, lvls, lvls, 1);
          fprintf(OutPutFile,'L1: %10.4e / %10.4e \n', norm(Afun(f)-RSS_apply_lu(G_A,f))/norm(Afun(f)), norm(Afun(f)-hodlr_apply(HODLR_A,f))/norm(Afun(f)))
          [Y_A, YB_A, YC_A, T_A, R_A, rk_A] = hodlrqr(HODLR_A, [], [], [], lvls, 1, tol_RSS);
          b = apply_bf(Factor, f);
        end
        for tol=[1E-10]

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
              if strcmp(initf, 'sol')
                f0 = RSS_inv_lu(G, N*ifft(apply_bf_adj(Factor, b)));
              elseif strcmp(initf, 'randn')
                f0 = randn(N,1);
              elseif strcmp(initf, 'rand')
                f0 = rand(N,1);
              end
              [x, iter] = SplitBregman(regm, mu, lambda, 1/N, @(x)N*ifft(apply_bf_adj(Factor,x))/N, @(x)hodlr_apply(HODLR_A, x), b, @(x)gradientxp_rc(x, 1/N)/N, @(x)gradientxp_rc(x, 1/N)/N, maxit1, maxit2, tol, @(x)RSS_inv_lu(G_A, x), f0, f0, zeros(N,1));
            elseif strcmp(regm,'F-L1')
              if strcmp(initf, 'sol')
                f0 = RSS_inv_lu(G, N*ifft(apply_bf_adj(Factor, b)));
              elseif strcmp(initf, 'randn')
                f0 = randn(N,1);
              elseif strcmp(initf, 'rand')
                f0 = rand(N,1);
              end
              [x, iter] = SplitBregman(regm, mu, lambda, 1/N, @(x)N*ifft(apply_bf_adj(Factor,x)), @(x)hodlr_apply(HODLR_A, x), b, @(x)N*ifft(x), @(x)fft(x), maxit1, maxit2, tol, @(x)RSS_inv_lu(G_A, x), f0, f0, zeros(N,1));
            elseif strcmp(regm,'L1')
              if strcmp(initf, 'sol')
                f0 = RSS_inv_lu(G, apply_bf_adj(Factor, b));
              elseif strcmp(initf, 'randn')
                f0 = randn(N,1);
              elseif strcmp(initf, 'rand')
                f0 = rand(N,1);
              end
              [x, iter] = SplitBregman(regm, mu, lambda, 1/N, @(x)apply_bf_adj(Factor,x), @(x)hodlr_apply(HODLR_A, x), b, @(x)x, @(x)x, maxit1, maxit2, tol, @(x)RSS_inv_lu(G_A,x), f0, f0, zeros(N,1));
            end
            t = toc;
            relerr = norm(f-x)/norm(f);
            fprintf(fileID,'Split Bregman with HIF guess           : tol %10.2e,   time/#iter1/#iter2 %10.4e / %4d / %4d, relerr %10.4e \n',tol,t,iter(1), iter(2), relerr);
            fprintf(OutPutFile, 'Split Bregman with HIF guess      : tol %10.2e,   time/#iter1/#iter2 %10.4e / %4d / %4d, relerr %10.4e \n',tol,t,iter(1), iter(2), relerr);
            solerrpcg_hif(i) = relerr;
            iters_hif(i) = iter(1)*iter(2);
            

            tic
            % [x,flag,relres,iter] = pcg(@(x) apply_bf_adj(Factor,apply_bf(Factor,x)),b,tol,maxit,@(x) hodlrqr_inv(Y, T, R, x)); 
            % [x, iter] = LinearizedBregman(@(x)apply_bf_adj(Factor, x), @(x)apply_bf(Factor, x), b, delta, mu, tol, maxit, apply_bf(Factor, hodlrqr_inv(Y, T, R, b)));
            % [x, iter] = SplitBregman(@(x)apply_bf_adj(Factor, x), @(x)apply_bf(Factor, x), b, @(x)x, @(x)eye(N), lambda, tol, maxit, apply_bf(Factor, hodlrqr_inv(Y, T, R, b)), apply_bf(Factor, hodlrqr_inv(Y, T, R, b)));
            if strcmp(regm,'TV-L1')
              if strcmp(initf, 'sol')
                f0 = hodlrqr_inv(Y, T, R, N*ifft(apply_bf_adj(Factor, b)));
              elseif strcmp(initf, 'randn')
                f0 = randn(N,1);
              elseif strcmp(initf, 'rand')
                f0 = rand(N,1);
              end
              [x, iter] = SplitBregman(regm, mu, lambda, 1/N, @(x)N*ifft(apply_bf_adj(Factor,x))/N, @(x)hodlr_apply(HODLR_A, x), b, @(x)gradientxp_rc(x, 1/N)/N, @(x)gradientxp_rc(x, 1/N)/N, maxit1, maxit2, tol, @(x)hodlrqr_inv(Y_A, T_A, R_A, x), f0, f0, zeros(N,1));
              % fprintf('TV-L1 cond: %10.4e / %10.4e \n', cond(hodlrqr_inv(Y_A, T_A, R_A, hodlr_apply(HODLR_A, eye(N)))), cond(hodlr_apply(HODLR_A, eye(N))))
            elseif strcmp(regm,'F-L1')
              if strcmp(initf, 'sol')
                f0 = hodlrqr_inv(Y, T, R, N*ifft(apply_bf_adj(Factor, b)));
              elseif strcmp(initf, 'randn')
                f0 = randn(N,1);
              elseif strcmp(initf, 'rand')
                f0 = rand(N,1);
              end
              [x, iter] = SplitBregman(regm, mu, lambda, 1/N, @(x)N*ifft(apply_bf_adj(Factor,x)), @(x)hodlr_apply(HODLR_A, x), b, @(x)N*ifft(x), @(x)fft(x), maxit1, maxit2, tol, @(x)hodlrqr_inv(Y_A, T_A, R_A, x), f0, f0, zeros(N,1));
              % fprintf('F-L1 cond: %10.4e / %10.4e \n', cond(hodlrqr_inv(Y_A, T_A, R_A, hodlr_apply(HODLR_A, eye(N)))), cond(hodlr_apply(HODLR_A, eye(N))))
            elseif strcmp(regm,'L1')
              if strcmp(initf, 'sol')
                f0 = hodlrqr_inv(Y, T, R, apply_bf_adj(Factor, b));
              elseif strcmp(initf, 'randn')
                f0 = randn(N,1);
              elseif strcmp(initf, 'rand')
                f0 = rand(N,1);
              end
              [x, iter] = SplitBregman(regm, mu, lambda, 1/N, @(x)apply_bf_adj(Factor,x), @(x)hodlr_apply(HODLR_A, x), b, @(x)x, @(x)x, maxit1, maxit2, tol, @(x)hodlrqr_inv(Y_A, T_A, R_A, x), f0, f0, zeros(N,1));
              % fprintf('L1 cond: %10.4e / %10.4e \n', cond(hodlrqr_inv(Y_A, T_A, R_A, hodlr_apply(HODLR_A, eye(N)))), cond(hodlr_apply(HODLR_A, eye(N))))
            end
            t = toc;
            relerr = norm(f-x)/norm(f);
            fprintf(fileID,'Split Bregman with HQR guess           : tol %10.2e,   time/#iter1/#iter2 %10.4e / %4d / %4d, relerr %10.4e \n',tol,t,iter(1), iter(2), relerr);
            fprintf(OutPutFile, 'Split Bregman with HQR guess      : tol %10.2e,   time/#iter1/#iter2 %10.4e / %4d / %4d, relerr %10.4e \n',tol,t,iter(1), iter(2), relerr);
            solerrpcg_hqr(i) = relerr;
            iters_hqr(i) = iter(1)*iter(2);
        end
      end
    end

end

