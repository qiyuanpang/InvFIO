function [FactorT,Factor, relerr, ZL, ZU] = run_bf_explicit_mod(N, Afun, mR, tol, fid)


    data_path = './data/';
    
    if(~exist(data_path, 'dir'))
        mkdir(data_path);
    end
    
    f = randn(N,1) + 1i*randn(N,1);
    
    nps = max(3,ceil(log2(mR)));
    % nps = mR;
    lsz = (2^nps)^2;

    
    tic;
    % Factor = bf_explicit(fun, xx, xbox, kk, kbox, mR, tol, 1);
    [Factor,ZL,ZU] = HSSBF_RS_fwd(Afun,(1:N)',(1:N)',mR,tol,lsz,1);
    FactorT = toc;
    
    tic;
    yy = HSSBF_apply(Factor,f);
    % yy = apply_bf(Factor, f);
    ApplyT = toc;
    RunT = FactorT + ApplyT;
    
    % norm(yy-fun(xx,kk)*f)/norm(fun(xx,kk)*f)
    
    NC = 256;
    tic;
    relerr = bf_explicit_check(N,Afun,f,(1:N)',(1:N)',yy,NC);
    Td = toc;
    Td = Td*N/NC;
    
    
    fprintf(fid,'------------------------------------------\n');
    fprintf(fid,'N                 : %4d\n', N);
    fprintf(fid,'Max Rank          : %4d\n', mR);
    fprintf(fid,'Tolerance         : %.3e\n', tol);
    fprintf(fid,'Relative Error_2  : %.3e\n', relerr);
    fprintf(fid,'Direct Time       : %.3e s\n', Td);
    fprintf(fid,'Running Time      : %.3e mins\n', RunT/60);
    fprintf(fid,'Factorization Time: %.3e mins\n', FactorT/60);
    fprintf(fid,'Applying Time     : %.3e s\n', ApplyT);
    fprintf(fid,'------------------------------------------\n\n');
    
    %[x,flag,relres,iter] = gmres(@(f) apply_bf(Factor,f),rand(N,1),[],10^(-12),N);
    %fprintf(fid,'Number iterations GMRES: %4d\n', iter(2));
    %tic;
    %[x,flag,relres,iter] = gmres(@(f) apply_bf(Factor,f),rand(N,1),[],10^(-12),N,@(f) apply_bf_adj(Factor,f));
    %GMRES_At = toc;
    %fprintf(fid,'Number GMRES iterations and time with A* preconditioning : %4d, %.3e s\n', iter(2),GMRES_At);
    
    %save([data_path 'Factor_' func_name '_' num2str(N) '_' num2str(mR) '.mat'],'Factor','-v7.3');
    
    end
