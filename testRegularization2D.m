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

initf = 'randn'
regm = 'F-L1'
whichlambda = 'lmus'

fun = @fun2;
func_name = 'fun2';
OutPutFile = fopen(['comp_2d/Regularization2D_save1_',func_name,'_',regm,'_',initf,'_',whichlambda,'.txt'],'w');


mR = 20;
%maxRank=[10 10 10];
mRk = 6;
%occ = 4*4;
tol_bf = 1E-12;
tol_peel = 1E-12;
tol_HIF = 1E-12;
maxit1 = 40;
maxit2 = 20;
repeat_num = 1;
delta = 10;
%mus = [2^(-12) 2^(-10) 2^(-8) 2^(-6) 2^(-4) 2^(-2) 2^0 2^2 2^4 2^6 2^8 2^(10) 2^(12)];
mus = [2^8];
lambdas = [2^(-16) 2^(-12) 2^(-8) 2^(-4) 2^0 2^4 2^(8) 2^(12) 2^(16)];



%dims = 2.^[14 15 16 17 18]
dims = [32 64 128 256 512]
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
    kN = log2(N);
    if kN < 6
       occ = 4*4;
    else
       occ = 8*8;
    end
    nlvl = floor(log2(N/sqrt(occ))) + 1;
    maxRank = floor(ones(1)*mRk*kN);
    fileID = fopen(['results_2d/MBF_',num2str(N),'.txt'],'w');

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

    fprintf(fileID,'\n');
    fprintf(fileID,'\n');
    fprintf(fileID,'------------------------------------------\n');
    fprintf(fileID,'N: %10d \n',N);
    
    fprintf(OutPutFile, '\n');
    fprintf(OutPutFile, '\n');
    fprintf(OutPutFile, '------------------------------------------\n');
    fprintf(OutPutFile, 'N: %10d \n',N);

    %% COnstruct BF factorization
    % [bftime(i), Factor, bferr(i)] = run_bf_explicit(N, func_name, mR, tol_bf, fileID);
    tic;
    Factor = mbf_explicit(fun, xx, xbox, kk, kbox, mR, tol_bf, 1);
    bftime(i) = toc;
    fprintf(OutPutFile, 'BF Fac err/time: %10.4e/%10.4e (s) \n', bferr(i), bftime(i));


    %% Construct HODLR factorization using peeling algorithm
    tStart=tic;
    for j = 1:repeat_num
      if strcmp(regm, 'L1')
        Afun = @(x) apply_mbf_adj(Factor,apply_mbf(Factor, x));
        [F, obranks, T, varargout] = peel_strong(Afun, kk', occ, tol_peel,maxRank,fileID);
      elseif strcmp(regm, 'TV-L1')
        Afun = @(x)  reshape(N^2*ifft2(reshape(apply_mbf_adj(Factor,apply_mbf(Factor, reshape(fft2(reshape(x,N,N,[])),N^2,[]))), N,N,[])),N^2,[]);
        [F, obranks, T, varargout] = peel_strong(Afun, kk', occ, tol_peel,maxRank,fileID);
      elseif strcmp(regm, 'F-L1')
        % Afun = @(x) N*ifft(apply_mbf_adj(Factor,apply_mbf(Factor, fft(x))));
        Afun = @(x)  reshape(N^2*ifft2(reshape(apply_mbf_adj(Factor,apply_mbf(Factor, reshape(fft2(reshape(x,N,N,[])),N^2,[]))), N,N,[])),N^2,[]);
        [F, obranks, T, varargout] = peel_strong(Afun, kk', occ, tol_peel,maxRank,fileID);
      end
    end
    t = toc(tStart)/repeat_num;
    factime_hodlr(i) = t;
    fprintf(fileID,'time peeling_strong: %10.4e \n',t);
    fprintf(OutPutFile,'time peeling_strong: %10.4e \n',t);
    fprintf(OutPutFile, 'observedMaxRanks: ');
    for ii = 1:length(varargout)
        fprintf(OutPutFile, '%6d / ', varargout(ii))
    end
    fprintf(OutPutFile, '\n');
    
    if strcmp(regm, 'L1')
      f = zeros(N^2,1);
      f(N^2/4) = 1;
    elseif strcmp(regm, 'TV-L1')
      f = ones(N^2,1);
      f(N^2/2+1:N^2) = 2;
    elseif strcmp(regm, 'F-L1')
      f1 = zeros(N^2,1);
      f1(N^2/4) = 1;
      f = ifft(f1);
    end
    % sum(abs(f) ~= 0)/N
    
    % Construct RSS factorization of HODLR matrix
    opts = struct('skip',1,'symm','h','verb',1);
    tStart_HIF = tic;
    for j = 1:repeat_num
    %   [G] = RSS_ldl(F,tol_RSS,fileID);
        [G,T2] = HIF_Hmatrix(F,kk',tol_HIF,T,opts,fileID);
    end
    t = toc(tStart_HIF)/repeat_num;
    factime_hif(i) = t;
    fprintf(fileID,'time HIF construction %10.4e \n',t);
    fprintf(OutPutFile, 'HIF Fac time: %10.4e \n', t);
   
    tic;
    for j = 1:repeat_num
      hifie_mv(G,f);
    end
    t=toc/repeat_num;
    % [e,niter] = snorm(N,@(x)(apply_bf_adj(Factor,apply_bf(Factor,x)) - RSS_apply(G,x)),[],[],32);
    % e = e/snorm(N,@(x)(apply_bf_adj(Factor,apply_bf(Factor,x))),[],[],1);
    if strcmp(regm, 'L1')
      e = norm(Afun(f) - hifie_mv(G,f))/norm(Afun(f));
    elseif strcmp(regm, 'TV-L1')
      e = norm(Afun(f) - hifie_mv(G,f))/norm(Afun(f));
    elseif strcmp(regm, 'F-L1')
      e = norm(Afun(f) - hifie_mv(G,f))/norm(Afun(f));
    end
    fprintf(fileID,'mv: %10.4e  time %10.4e\n',e,t);
    apptime_hif(i) = t;
    apperr_hif(i) = e;
    fprintf(OutPutFile, 'HIF mv err/time: %10.4e/%10.4e \n', e, t);

    tic;
    for j = 1:repeat_num
      hifie_sv(G,f);
    end
    t=toc;
    % NORM(INV(A) - INV(F))/NORM(INV(A)) <= NORM(I - A*INV(F))
    %[e,niter] = snorm(N,@(x)(x-RSS_inv(G,apply_bf_adj(Factor,apply_bf(Factor,x)))),[],[],1);
    %fprintf(fileID,'sv: %10.4e / (%4d) time %10.4e\n',e,niter,t);
    %NORM(INV(L*L') - INV(F))/NORM(INV(L*L')) <= NORM(I - L*INV(F)*L')
    % [e,niter] = snorm(N,@(x)(x-apply_bf(Factor,RSS_inv(G,apply_bf_adj(Factor,x)))),[],[],32);
    if strcmp(regm, 'L1')
      e = norm(f-Afun(hifie_sv(G, f)))/norm(f);
    elseif strcmp(regm, 'TV-L1')
      e = norm(f-Afun(hifie_sv(G, f)))/norm(f);
    elseif strcmp(regm, 'F-L1')
      e = norm(f-Afun(hifie_sv(G, f)))/norm(f);
    end
    fprintf(fileID,'sv: %10.4e time %10.4e\n',e,t);
    soltime_hif(i) = t;
    solerr_hif(i) = e;
    fprintf(OutPutFile, 'HIF sv err/time: %10.4e/%10.4e \n', e, t);

    % run SB
    fprintf(OutPutFile, 'regularizer: %s  init-scheme: %s \n', regm, initf)
    for mu = mus
      %lmus = [mu/128 mu/64 mu/32 mu/16 mu/8 mu/4 mu/2 mu mu*2 mu*4 mu*8 mu*16 mu*32 mu*64 mu*128];
      lmus = [mu/32];
      if strcmp(whichlambda,'ind')
        LAMBDAS = lambdas;
      elseif strcmp(whichlambda, 'lmus');
        LAMBDAS = lmus;
      end
      for lambda = LAMBDAS
        fprintf(OutPutFile, 'mu: 2^{%4d}  lambda: 2^{%4d} \n', log2(mu), log2(lambda))
        if strcmp(regm,'TV-L1')
          Afun = @(x)  reshape(mu*N^2*ifft2(reshape(apply_mbf_adj(Factor,apply_mbf(Factor, reshape(fft2(reshape(x,N,N,[])),N^2,[]))), N,N,[])),N^2,[])+...
                 lambda*reshape(laplacianp_rc(reshape(x,N,N,[]),1/N),N^2,[])+lambda*reshape(pagetranspose(laplacianp_rc(pagetranspose(reshape(x,N,N,[])),1/N)),N^2,[]);
        %   Afun = @(x) mu*N*ifft(apply_mbf_adj(Factor,apply_mbf(Factor,fft(x)))) + lambda*laplacianp_rc(x,1/N);
          fprintf(fileID, 'Information for TV-L1 \n')
          [F_A, obranks, T, varargout] = peel_strong(Afun, kk', occ, tol_peel,maxRank,fileID);
          [G_A,T2] = HIF_Hmatrix(F_A,kk',tol_HIF,T,opts,fileID);
          fprintf(OutPutFile, 'SB TV-L1: %10.4e / %10.4e \n', norm(Afun(f)-hifie_mv(G_A,f))/norm(Afun(f)), norm(f-Afun(hifie_sv(G_A,f)))/norm(f))
          b = apply_mbf(Factor,reshape(fft2(reshape(f,N,N,[])),N^2,[]));
        elseif strcmp(regm,'F-L1')
          Afun = @(x)  reshape(mu*N^2*ifft2(reshape(apply_mbf_adj(Factor,apply_mbf(Factor, reshape(fft2(reshape(x,N,N,[])),N^2,[]))), N,N,[])),N^2,[])+lambda*N^2*x;
          fprintf(fileID, 'Information for F-L1 \n')
          [F_A, obranks, T, varargout] = peel_strong(Afun, kk', occ, tol_peel,maxRank,fileID);
          [G_A,T2] = HIF_Hmatrix(F_A,kk',tol_HIF,T,opts,fileID);
          fprintf(OutPutFile, 'SB F-L1: %10.4e / %10.4e \n', norm(Afun(f)-hifie_mv(G_A,f))/norm(Afun(f)), norm(f-Afun(hifie_sv(G_A,f)))/norm(f))
          b = apply_mbf(Factor,reshape(fft2(reshape(f,N,N,[])),N^2,[]));
        elseif  strcmp(regm,'L1')
          Afun = @(x) mu*apply_mbf_adj(Factor, apply_mbf(Factor, x)) + lambda*x;
          fprintf(fileID, 'Information for L1 \n')
          [F_A, obranks, T, varargout] = peel_strong(Afun, kk', occ, tol_peel,maxRank,fileID);
          [G_A,T2] = HIF_Hmatrix(F_A,kk',tol_HIF,T,opts,fileID);
          fprintf(OutPutFile, 'SB L1: %10.4e / %10.4e \n', norm(Afun(f)-hifie_mv(G_A,f))/norm(Afun(f)), norm(f-Afun(hifie_sv(G_A,f)))/norm(f))
          b = apply_mbf(Factor, f);
        end
        fprintf(OutPutFile, 'observedMaxRanks: ');
        for ii = 1:length(varargout)
            fprintf(OutPutFile, '%6d / ', varargout(ii))
        end
        fprintf(OutPutFile, '\n');

        for tol=[1E-10]

            tic
            if strcmp(regm,'TV-L1')
              if strcmp(initf, 'sol')
                f0 = hifie_sv(G, reshape(N^2*ifft2(reshape(apply_mbf_adj(Factor, b), N,N,[])),N^2,[]));
              elseif strcmp(initf, 'randn')
                f0 = randn(N^2,1);
              elseif strcmp(initf, 'rand')
                f0 = rand(N^2,1);
              end
              norm(f0-f)/norm(f)
              [x, iter] = SplitBregman_TV2D(mu, lambda, @(x)reshape(N^2*ifft2(reshape(apply_mbf_adj(Factor, b), N,N,[])),N^2,[]), @(x)hifie_mv(G_A,x), b, @(x)reshape(gradientxp_rc(reshape(x,N,N,[]), 1/N),N^2,[]),...
                          @(x)reshape(pagetranspose(gradientxp_rc(pagetranspose(reshape(x,N,N,[])), 1/N)),N^2,[]), @(x)reshape(gradientxp_rc(reshape(x,N,N,[]), 1/N),N^2,[]), @(x)reshape(pagetranspose(gradientxp_rc(pagetranspose(reshape(x,N,N,[])), 1/N)),N^2,[]), maxit1, maxit2, tol, @(x)hifie_sv(G_A, x), f0, f0, zeros(N^2,1),f0,zeros(N^2,1));
            elseif strcmp(regm,'F-L1')
              if strcmp(initf, 'sol')
                f0 = hifie_sv(G, reshape(N^2*ifft2(reshape(apply_mbf_adj(Factor, b), N,N,[])),N^2,[]));
              elseif strcmp(initf, 'randn')
                f0 = randn(N^2,1);
              elseif strcmp(initf, 'rand')
                f0 = rand(N^2,1);
              end
            %   norm(f0-f)/norm(f)
              [x, iter] = SplitBregman(regm, mu, lambda, 1/N, @(x)reshape(N^2*ifft2(reshape(apply_mbf_adj(Factor, b), N,N,[])),N^2,[]), @(x)hifie_mv(G_A, x), b, @(x)reshape(N^2*ifft2(reshape(x,N,N,[])),N^2,[]), @(x)reshape(fft2(reshape(x,N,N,[])),N^2,[]), maxit1, maxit2, tol, @(x)hifie_sv(G_A, x), f0, f0, zeros(N^2,1));
            elseif strcmp(regm,'L1')
              if strcmp(initf, 'sol')
                f0 = hifie_sv(G, apply_mbf_adj(Factor, b));
              elseif strcmp(initf, 'randn')
                f0 = randn(N^2,1);
              elseif strcmp(initf, 'rand')
                f0 = rand(N^2,1);
              end
              [x, iter] = SplitBregman(regm, mu, lambda, 1/N, @(x)apply_mbf_adj(Factor,x), @(x)hifie_mv(G_A, x), b, @(x)x, @(x)x, maxit1, maxit2, tol, @(x)hifie_sv(G_A,x), f0, f0, zeros(N^2,1));
            end
            t = toc;
            relerr = norm(f-x)/norm(f);
            fprintf(fileID,'Split Bregman with HIF guess           : tol %10.2e,   time/#iter1/#iter2  %10.4e / %4d / %4d, relerr %10.4e \n',tol,t,iter(1), iter(2), relerr);
            fprintf(OutPutFile, 'Split Bregman with HIF guess      : tol %10.2e,   time/#iter1/#iter2  %10.4e / %4d / %4d, relerr %10.4e \n',tol,t,iter(1), iter(2), relerr);
            solerrpcg_hif(i) = relerr;
            iters_hif(i) = iter(1)*iter(2);
            
        end
      end
    end

end
