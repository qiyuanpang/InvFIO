clear all
startup; 
curpath = pwd;
dirs = {'BF/1D'};
for s = dirs
  addpath(sprintf('%s/%s/src',curpath,s{:}));
  addpath(sprintf('%s/%s/test',curpath,s{:}));
end


func_name = 'fun_FIO_var1';%'fun_FIO_var2';'fun_FIO';
  

mR = 8;
occ = 32;
tol_bf = 10^(-6);
tol_peel = 10^(-4);
tol_RSS = 10^(-3);

N_size = [];
t_BFF =[];
t_HODLR = [];
t_HIF = [];
t_apply = [];

for N= [256,1024,4096,2*8192,8*8192,32*8192]

    N_size = [N_size,N];
    
    fileID = fopen(['results_1d/BFF_',func_name,'_N_',num2str(N),'.txt'],'w');

    fprintf(fileID,'\n');
    fprintf(fileID,'\n');
    fprintf(fileID,'------------------------------------------\n');
    fprintf(fileID,'N: %10d \n',N);
    
    %Factorize FIO matrix with BFF
    [t,Factor] = run_bf_explicit(N, func_name, mR, tol_bf, fileID);
    t_BFF =[t_BFF,t];
    
    %Construct HODLR matrix
    tStart=tic;
    [F,HOLDR] = HODLR_construction( N, @(x) apply_bf(Factor,x), @(x) apply_bf_adj(Factor,x), tol_peel, fileID, occ, 20,10);
    t = toc(tStart);
    t_HODLR = [t_HODLR,t];
    fprintf(fileID,'time H_matrix construction 10e-6: %10.4e )\n',t);
    % NORM(A - F)/NORM(A)
    [e,niter] = snorm2(N,@(x)(apply_bf_adj(Factor,apply_bf(Factor,x)) - HODLR_apply(HOLDR,x)),[],[],1);
    e = e/snorm2(N,@(x)(apply_bf_adj(Factor,apply_bf(Factor,x))),[],[],1);
    fprintf(fileID,'mv: %10.4e / %4d )\n',e,niter);
    
    
    % Construct RSS factorization of HODLR matrix
    tStart_HIF = tic;
    [G] = RSS(F,tol_RSS,fileID);
    t = toc(tStart_HIF);
    t_HIF = [t_HIF,t];
    fprintf(fileID,'time HIF construction %10.4e \n',t);
    
    f = randn(N,1) + 1i*randn(N,1);
    tic
    RSS_apply(G,f);
    t=toc;
    [e,niter] = snorm2(N,@(x)(apply_bf_adj(Factor,apply_bf(Factor,x)) - RSS_apply(G,x)),[],[],1);
    e = e/snorm2(N,@(x)(apply_bf_adj(Factor,apply_bf(Factor,x))),[],[],1);
    fprintf(fileID,'mv: %10.4e / (%4d ) time %10.4e\n',e,niter,t);
    tic
    RSS_inv(G,f);
    t=toc;
    % NORM(INV(A) - INV(F))/NORM(INV(A)) <= NORM(I - A*INV(F))
    %[e,niter] = snorm2(N,@(x)(x-RSS_inv(G,apply_bf_adj(Factor,apply_bf(Factor,x)))),[],[],1);
    %fprintf(fileID,'sv: %10.4e / (%4d) time %10.4e\n',e,niter,t);
    %NORM(INV(L*L') - INV(F))/NORM(INV(L*L')) <= NORM(I - L*INV(F)*L')
    [e,niter] = snorm2(N,@(x)(x-apply_bf(Factor,RSS_inv(G,apply_bf_adj(Factor,x)))),[],[],1);
    fprintf(fileID,'sv: %10.4e / (%4d) time %10.4e\n',e,niter,t);
    
    
    tic
    RSS_inv(G,apply_bf_adj(Factor,f));
    t =toc;
    t_apply = [t_apply,t];
   
    fprintf(fileID,'inv(AtA)At time: %10.4e\n',t);
    tic
    apply_bf_adj(Factor,f);
    t = toc;
    fprintf(fileID,'At time: %10.4e\n',t);
    
    % run CG 
     for tol=[10^(-4),10^(-8)]
         
         tic
         [~,flag,relres,iter] = pcg(@(x) apply_bf_adj(Factor,apply_bf(Factor,x)),apply_bf_adj(Factor,apply_bf(Factor,f)),tol,N);
         t1 = toc;
         fprintf(fileID,'CG with A* preconditioning: tol %10.2e , time/#iter %10.4e / %4d )\n',tol,t1,iter);
         tic
         [~,flag,relres,iter] = pcg(@(x) apply_bf_adj(Factor,apply_bf(Factor,x)),apply_bf_adj(Factor,apply_bf(Factor,f)),tol,N,@(x) RSS_inv(G,x)); 
          t = toc;
         fprintf(fileID,'CG with HIF preconditioning: tol %10.2e,   time/#iter %10.4e / %4d )\n',tol,t,iter);
     end
     
end

