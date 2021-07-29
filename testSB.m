clear all;
dims = [13];
regm = 'L1';
mus = [2^(-10) 2^(-6) 2^1 2^6 2^(10)];
lambdas = [1E-12 1E-9 1E-4 1E0 1E4 1E9 1E12];
tol = 1E-15;
maxit1 = 10;
maxit2 = 50;
OutPutFile = fopen(['testSB.txt'],'w');
soltype = 3;
for i = dims
    N = 2^i;
    M = 10;
    F = randn(M,N);
    if soltype == 1
        f = zeros(N,1);
        f(N/4) = 1;
        b = fft(f);
        b = b(1:N/2);
    elseif soltype == 2        
        f = ones(N,1);
        f(N/4+1:N/2) = 2;
        b = fft(f);
        b = b(1:N/2);
    elseif soltype == 3
        f = zeros(N,1);
        f(N/4) = 1;
        b = F*f;
        % b = b(1:M);
    end
    hf = @(x)x(1:N/2);
    hfM = @(x)x(1:M);
    ft = [ones(N/2,1);zeros(N/2,1)];
    ftM = [ones(M,1);zeros(N-M,1)];
    f0s = [f 2*f randn(N,1)];
    for mu = mus
        for lambda = lambdas
            data = [];
            for kk = 1:size(f0s,2)
                f0 = f0s(:,kk);
                % [x, iter] = SplitBregman('L1', mu, lambda, 1/N, @(x)N*ifft(x), @(x)mu*N*ifft(fft(x))+lambda*x, b, @(x)x, @(x)x, maxit1, maxit2, tol, @(x)ifft(fft(x))/(N*mu+lambda), f0, f0, zeros(N,1));
                % [x, iter] = SplitBregman('L1', mu, lambda, 1/N, @(x)N*ifft([x;zeros(N/2,1)]), @(x)mu*N*ifft(fft(x).*ft)+lambda*(x.*ft), b, @(x)[x;zeros(N/2,1)], @(x)x(1:N/2), maxit1, maxit2, tol, @(x)x, f0, f0(1:N/2), zeros(N/2,1));
                [x, iter] = SplitBregman('L1', mu, lambda, 1/N, @(x)F'*x, @(x)mu*F'*F*x+lambda*x.*ftM, b, @(x)[x;zeros(N-M,1)], @(x)x(1:M), maxit1, maxit2, tol, @(x)x, f0, f0(1:M), zeros(M,1));
                % [x, iter] = SplitBregman('TV-L1', mu, lambda, 1/N, @(x)N*ifft(x), @(x)mu*N*ifft(fft(x))+lambda*laplacianp(x, 1/N^2), b, @(x)gradientxp(x,1/N), @(x)x, maxit1, maxit2, tol, @(x)ifft(fft(x))/mu/N, f0, f0,zeros(N,1));
                % [x, iter] = SplitBregman('TV-L1', mu, lambda, 1/N, @(x)N*ifft([x;zeros(N/2,1)]), @(x)mu*N*ifft(fft(x).*ft)+lambda*laplacianp(x, 1/N^2).*ft, b, @(x)[gradientxp(x,1/N);zeros(N/2,1)], @(x)hf(gradientxp(x,1/N)), maxit1, maxit2, tol, @(x)x, f0, f0(1:N/2),zeros(N/2,1));
                % [x, iter] = SplitBregman('TV-L1', mu, lambda, 1/N, @(x)F'*x, @(x)mu*F'*F*x+lambda*laplacianp(x,1/N^2).*ftM, b, @(x)[gradientxp(x,1/N);zeros(N-M,1)], @(x)hfM(gradientxp(x,1/N)), maxit1, maxit2, tol, @(x)x, f0, f0(1:M),zeros(M,1));
                
                % [x, iter] = SplitBregman('TV-L1', mu, lambda, 1/N, @(x)N*ifft(x), @(x)mu*N*ifft(fft(x))+lambda*laplacianp(x, 1/N^2).*gradientxp(x,1/N), b, @(x)laplacianp(x,1/N^2), maxit1, maxit2, tol, @(x)ifft(fft(x))/mu/N, f, f,zeros(N,1));
                % [x, iter] = SplitBregman('TV-L1', mu, lambda, 1/N, @(x)N*ifft(x), @(x)mu*N*ifft(fft(x))+lambda*sum(laplacianp(x, 1/N^2)), b, @(x)gradientxpT(x,1/N), maxit1, maxit2, tol, @(x)ifft(fft(x))/mu/N, f, f,zeros(N,1));
                % [x, iter] = SplitBregman('TV-L1', mu, lambda, 1/N, @(x)N*ifft(x), @(x)mu*N*ifft(fft(x))+lambda*gradientxpT(gradientxp(x,1/N),1/N), b, @(x)gradientxpT(x,1/N), maxit1, maxit2, tol, @(x)ifft(fft(x))/mu/N, f, f,zeros(N,1));
                if soltype < 3
                    y = fft(x);
                    e = norm(b-y(1:N/2))/norm(b);
                else
                    y = F*x;
                    e = norm(b-y)/norm(b);
                end
                % e = norm(f-x)/norm(f);
                data = [data iter(1) iter(2) e];
            end
            fprintf(OutPutFile,'\\hline \n')
            fprintf(OutPutFile,'($2^{%d }$, $10^{%d }$) & %d & %d & %5.1e & %d & %d & %5.1e & %d & %d & %5.1e \\\\ \n', log2(mu), log10(lambda), data(1), data(2), data(3), data(4), data(5), data(6), data(7), data(8), data(9))
        end
    end
    % fprintf('results %10.4e / %d / %d \n', e, iter(1), iter(2))

end