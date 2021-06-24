close all;
clear all;

startup;

tol = 1e-12
kernel = "7"
fun_name = "fun" + kernel;
func_name = fun_name;
occ = 32;
rank_or_tol = 1e-3
tol_sol = 1e-8
maxit = 200
repeat_num = 5;
n0 = 8;
tt = 6;
rand_or_cheb = 'cheb';

% dims = [16 25 36]
dims = [16 25 36 49 64 81 100]
cases = length(dims);
apptime = zeros(cases, 1);
bferr = zeros(cases, 1);
condAs = zeros(cases, 1);
condATAs = zeros(cases, 1);

soltime = zeros(cases, 1);
apperr = zeros(cases, 1);
solerr = zeros(cases, 1);

ranks_hqr = zeros(cases, 1);
ranks_hodlr = zeros(cases, 1);
solerrpcg = zeros(cases, 1);
iters = zeros(cases, 1);
for i = 1:cases
    n = dims(i);
    N = n^2;
    NG = 5*floor(log2(N));
    rk = 9*floor(log2(N));
    k = -N/2:N/2-1;
    kk = k(:);
    %rank_or_tol = 8*floor(log2(N))    

    x = (0:N-1)/N;
    xx = x(:);
    

    fprintf('N = %4d \n', N)

    switch fun_name
        case "fun1"
            fun = @(x,k)fun_fio_1D(x,k);
        case "fun2"
            fun = @(x,k)fun_fio2_1D(x,k);
        case "fun3"
            fun = @(x,k)fun_fio3_1D(x,k,0.1);
        case "fun4"
            fun = @(x,k)fun_fio4_1D(N, x, k, 0.05);
        case "fun5"
            fun = @(x,k)fun_fio5_1D(x,k);
        case "fun6"
            fun = @(x,k)fun_fio6_1D(N, x, k, 0.1);
        case "fun7"
            fun = @(x,k)fun_fio6_1D(N, x, k, 0.05);

    end

    Factor = BF_IDBF(fun, xx, kk, n0, rk, tol, rand_or_cheb, tt);

    A = BF_apply(Factor,eye(N));
    ATA = BF_adj_apply(Factor, A);

    condA = cond(A);
    condATA = cond(ATA);
    fprintf('condition number estimation of A  : %10.4e \n', condA)
    fprintf('condition number estimation of ATA: %10.4e \n', condATA)

    condAs(i) = condA;
    condATAs(i) = condATA;
    
    A1 = fun(xx, kk);
    ATA1 = A1'*A1;

    condA1 = cond(A1);
    condATA1 = cond(ATA1);
    fprintf('condition number estimation of A1  : %10.4e \n', condA1)
    fprintf('condition number estimation of ATA1: %10.4e \n', condATA1)
    
    randx = rand(N, 1);
    bfacc = norm(A*randx-A1*randx)/norm(A1*randx);
    fprintf('err between matrices: %10.4e \n', bfacc)
    bferr(i) = bfacc;    


    lvls = round(log2(N/16));
    F = hodlr(ATA, lvls, 1E-9);
    % F
    % while ~isempty(F)
    %     F = F.A11
    % end
    ranks_hodlr(i) = F.maxrk;
    % norm2 = norm(ATA);
    norm2 = 1;
    tic;
    for j =1:repeat_num
        [Y, YB, YC, T, R, rk] = hodlrqr(F, [], [], [], lvls, norm2, rank_or_tol);
    end
    factime(i) = toc/repeat_num;
    ranks_hqr(i) = rk;

    
    v = rand(N, 1);
    tic;
    for j = 1:repeat_num
        y = hodlr_apply(R, v);
        y = y - hodlr_apply(Y, hodlr_apply(T, hodlr_adj_apply(Y, y)));
    end
    apptime(i) = toc/repeat_num;
    apperr(i) = norm(y-hodlr_apply(F, v))/norm(hodlr_apply(F, v));
    fprintf('hqr_mv   err: %10.4e/%10.4e (s) \n', apperr(i), apptime(i))

    tic;
    for j = 1:repeat_num
        x = y - hodlr_apply(Y, hodlr_adj_apply(T, hodlr_adj_apply(Y, y)));
        x = hodlr_tri_sol(R, x);
    end
    soltime(i) = toc/repeat_num;
    solerr(i) = norm(x-v)/norm(v);
    fprintf('hqr_sv   err/time: %10.4e/%10.4e (s) \n', solerr(i), soltime(i))
    
    pcd = @(y)hodlr_tri_sol(R, y - hodlr_apply(Y, hodlr_adj_apply(T, hodlr_adj_apply(Y, y))));
    % fun = @(v)hodlr_apply(R, v) - hodlr_apply(Y, hodlr_apply(T, hodlr_adj_apply(Y, hodlr_apply(R, v))));
    fun = @(v)hodlr_apply(F, v);
    y = fun(v);
    [x, flag, relres, iter] = pcg(fun, y, tol_sol, maxit, pcd);
    relerr = norm(x-v)/norm(v);
    solerrpcg(i) = relerr;
    iters(i) = iter;
    fprintf('Solve the equation by PCG with HQR as a preconditioner in %4d iterations, rel error: %10.4e \n', iter, relerr)
    
end

N = dims.^2;
logN = log2(N);
NlogN = logN + log2(logN);
N2logN = logN + 2*log2(logN);

fig = figure(1);
hold on;
h(1) = plot(logN, log2(apptime));
h(2) = plot(logN, NlogN-NlogN(1)+log2(apptime(1)));
h(3) = plot(logN, N2logN-N2logN(1)+log2(apptime(1)));
xlabel('Log(N)');
ylabel('Log(Time)/s'); 
title('Forward Time Scaling');
legend(h, 'App time', 'N log N', 'N log^2 N');
axis square;
saveas(fig, "./scaling/1D/kernel" + kernel + "/apptime_" + func_name + ".png");
hold off;

fig = figure(2);
hold on;
h(1) = plot(logN, log2(factime));
h(2) = plot(logN, NlogN-NlogN(1)+log2(factime(1)));
h(3) = plot(logN, N2logN-N2logN(1)+log2(factime(1)));
xlabel('Log(N)');
ylabel('Log(Time)/s'); 
title('Fac Time Scaling');
legend(h, 'Fac time', 'N log N', 'N log^2 N');
axis square;
saveas(fig, "./scaling/1D/kernel" + kernel + "/factime_" + func_name + ".png");
hold off;

fig = figure(3);
hold on;
h(1) = plot(logN, log10(apperr));
xlabel('Log(N)');
ylabel('Log10(error)'); 
title('Forward Error');
% legend(h, 'App time', 'N log N', 'N log^2 N');
axis square;
saveas(fig, "./scaling/1D/kernel" + kernel + "/apperr_" + func_name + ".png");
hold off;

fig = figure(4);
hold on;
h(1) = plot(logN, log10(solerr));
xlabel('Log(N)');
ylabel('Log10(error)'); 
title('Backward Error');
% legend(h, 'App time', 'N log N', 'N log^2 N');
axis square;
saveas(fig, "./scaling/1D/kernel" + kernel + "/solerr_" + func_name + ".png");
hold off;

fig = figure(5);
hold on;
h(1) = plot(logN, log2(ranks_hodlr));
h(2) = plot(logN, log2(ranks_hqr));
h(3) = plot(logN, 0.5*logN-0.5*logN(1)+(log2(ranks_hodlr(1))+log2(ranks_hqr(1)))/2);
h(4) = plot(logN, log2(logN)-log2(logN(1))+(log2(ranks_hodlr(1))+log2(ranks_hqr(1)))/2);
xlabel('Log(N)');
ylabel('Log2(rank)'); 
title('Rank');
legend(h, 'HODLR', 'HQR', 'sqrt(N)', 'log N');
saveas(fig, "./scaling/1D/kernel" + kernel + "/rank_"+ func_name + ".png");
hold off;

fig = figure(6);
hold on;
h(1) = plot(logN, log2(soltime));
h(2) = plot(logN, NlogN-NlogN(1)+log2(soltime(1)));
h(3) = plot(logN, N2logN-N2logN(1)+log2(soltime(1)));
xlabel('Log(N)');
ylabel('Log(Time)/s'); 
title('Backward Time Scaling');
legend(h, 'Backward time', 'N log N', 'N log^2 N');
axis square;
saveas(fig, "./scaling/1D/kernel" + kernel + "/soltime_" + func_name + ".png");
hold off;

fig = figure(7);
hold on;
h(1) = plot(logN, log10(solerrpcg));
xlabel('Log(N)');
ylabel('Log10(error)'); 
title('Error of PCG Solution');
% legend(h, 'HIF', 'HQR');
saveas(fig, "./scaling/1D/kernel" + kernel + "/pcgerr_"+ func_name + ".png");
hold off;

fig = figure(8);
hold on;
h(1) = plot(logN, iters);
% h(2) = plot(logN, iters_hqr);
xlabel('Log(N)');
ylabel('Iteration'); 
title('Iterations in PCG');
% legend(h, 'HIF', 'HQR');
saveas(fig, "./scaling/1D/kernel" + kernel + "/iters_" + func_name + ".png");
hold off;


exit;
