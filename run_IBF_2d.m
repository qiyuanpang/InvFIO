clear all
startup; 

curpath = pwd;
addpath(sprintf('%s/%s',curpath,'BF/2D/GBF/src'))
dirs = {'BF/2D/GBF'};
for s = dirs
  addpath(sprintf('%s/%s/src',curpath,s{:}))
  addpath(sprintf('%s/%s/test',curpath,s{:}))
end


fun = @fun0;
func_name = 'fun0';



occ = 8*8;
mR = 32;
tol_BFF = 10^(-6);
tol_peel = 10^(-4);
tol_HIF = 10^(-3);
maxRank = [18,18];

for kN=[1]
    
    N = 32*(2^kN);
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
    
    % Construct Butterfly factorization
    tic
    Factor = mbf_explicit(fun, xx, xbox, kk, kbox, mR, tol_BFF, 1);
    t = toc;
    fprintf(fileID,'BF time: %10.4e (s)\n',t);
    Afun = @(x) apply_mbf_adj(Factor,apply_mbf(Factor,x));
    
    % Construct Hierarchical matrix
    peel_time = tic; 
    [GHat, ranks, T, varargout] = peel_strong(Afun, kk', occ, tol_peel,[0,maxRank],fileID);
    t = toc(peel_time);
    fprintf(fileID,'Peeling time: %10.4e (s)\n',t);
    fprintf(fileID,['-'*ones(1,80) '\n']);
    
    
    % Invert Hierarchical matrix
    opts = struct('skip',1,'symm','h','verb',1);
    t_HIF = tic;
    [F,T2] = HIF_Hmatrix(GHat,kk',tol_HIF,T,opts,fileID);
    t = toc(t_HIF);
    fprintf(fileID,'HIF time: %10.4e (s)\n',t);
    
    
    f = randn(N^2,1) + 1i*randn(N^2,1);
    Bf = apply_mbf_adj(Factor,f);
    tic
    hifie_sv(F,Bf);
    t = toc;
    fprintf(fileID,'Apply HIF time: time %10.4e \n',t);
    
     % test accuracy using randomized power method
    X = rand(N^2,1);
    X = X/norm(X);
     % NORM(INV(A) - INV(F))/NORM(INV(A)) <= NORM(I - A*INV(F))
    tic
    Y = hifie_sv(F,X);
    t = toc;
    [e,niter] = snorm(N^2,@(x)(x - Afun(hifie_sv(F,x))),[],[],1);
    fprintf(fileID,'sv: %10.4e / %4d / %10.4e (s)\n',e,niter,t);
    %[e,niter] = snorm(N^2,@(x)(x - apply_mbf(Factor,hifie_sv(F,apply_mbf_adj(Factor,x)))),[],[],1);
    %fprintf(fileID,'inv A: %10.4e / %4d / %10.4e (s)\n',e,niter,t)



     %CG
     count = 0;
      for tol=[10^(-4),10^(-8)]
          count = count + 1;
          tic
          [x_1,flag,relres,iter] =  pcg(@(x) apply_mbf_adj(Factor,apply_mbf(Factor,x)),Bf,tol,N^2);
          t = toc;
          fprintf(fileID,'CG with A* preconditioning: tol %10.2e , time/#iter %10.4e / %4d )\n',tol,t,iter);%(2));

          tic
          [x_2,flag,relres,iter] =  pcg(@(x) (apply_mbf_adj(Factor,apply_mbf(Factor,x))),hifie_sv(F,Bf),tol,N^2,@(x) hifie_sv(F,x));
          t = toc;
          fprintf(fileID,'CG with HIF 0.001 preconditioning: tol %10.2e,   time/#iter %10.4e / %4d )\n',tol,t,iter); 
          
      end

end
