function [FactorT,Factor, relerr] = run_bf_explicit(N, func_name, mR, tol, fid)


data_path = './data/';

if(~exist(data_path, 'dir'))
    mkdir(data_path);
end

kbox = [1,N+1];
k = 1:N;
kk = k(:);
% noise = rand(N,1)*0.9;
% kk = kk + noise;

xbox = [1,N+1];
x = 1:N;
xx = x(:);
% noise = rand(N,1)*0.5;
% xx = xx + noise;
% xx(2) = xx(1) + 0.1;
% xx(N-1) = xx(N) - 0.1;
% xx(N/2+1) = xx(N/2) + 0.05;

switch func_name
    case 'fun_FIO'
        fun = @(x,k)fun0(N,x,k);
    case 'fun_FIO_var1'
        fun = @(x,k)fun0var1(N,x,k);
    case 'fun_FIO_var2'
        fun = @(x,k)fun0var2(N,x,k);
    case 'fun_FIO_period'
        fun = @(x,k)fun0period(N,x,k);
    case 'funF'
        fun = @(x,k)funF(N,x,k);
    case 'funH'
        fun = @(x,k)funH(N,x,k);
    case 'fun_FIO_5'
        fun = @(x,k)fun5(N,x,k);
    case 'fun_FIO_var4'
        fun = @(x,k)fun0var4(N,x,k);
    case 'fun_FIO_var6'
        fun = @(x,k)fun0var6(N,x,k);
    case 'fun_FIO_var7'
        fun = @(x,k)fun0var7(N,x,k);
    case 'fun_FIO_var8'
        fun = @(x,k)fun0var8(N,x,k);
    case 'fun_FIO_var9'
        fun = @(x,k)fun0var9(N,x,k);
    case 'fun_FIO_var10'
        fun = @(x,k)fun0var10(N,x,k);
    case 'fun_FIO_var11'
        fun = @(x,k)fun0var11(N,x,k);
    case 'fun_FIO_var12'
        fun = @(x,k)fun0var12(N,x,k);
end

f = randn(N,1) + 1i*randn(N,1);

nps = max(3,ceil(log2(mR)));
% nps = mR;
lsz = (2^nps)^2;

Afun = @(x,k)fun(xx(x), kk(k));

% fprintf('HSSBF error: %10.4e \n', norm(HSSBF_apply(Fac, f) - fun(xx,kk)*f)/norm(fun(xx,kk)*f));

% norm(HSSBF_apply(Fac,f) - HSSBF_apply(ZL, f) - HSSBF_apply(ZU,f)+f)
% norm(fun(xx,kk)*f-Afun((1:N)', (1:N)')*f)/norm(fun(xx,kk)*f)

tic;
Factor = bf_explicit(fun, xx, xbox, kk, kbox, mR, tol, 1);
% [Factor,ZL,ZU] = HSSBF_RS_fwd(Afun,(1:N)',(1:N)',mR,tol,lsz,1);
FactorT = toc;

tic;
% yy = HSSBF_apply(Factor,f);
yy = apply_bf(Factor, f);
ApplyT = toc;
RunT = FactorT + ApplyT;

% norm(yy-fun(xx,kk)*f)/norm(fun(xx,kk)*f)

NC = 256;
tic;
relerr = bf_explicit_check(N,fun,f,xx,kk,yy,NC);
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
