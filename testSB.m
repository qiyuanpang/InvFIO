clear all;
dims = [8];
regm = 'L1';
mus = [2^(-7) 2^(-4) 2^0 2^4 2^(7) ];
lambdas = [1E-12 1E-10 1E-8 1E-6 1E-4 1E-2 1E0 1E2 1E4];
% lambdas = [1E-2 1E2]
tol = 1E-13;
maxit1 = 50;
maxit2 = 50;
OutPutFile = fopen(['testSB.txt'],'w');
NoFIle = fopen(['less.txt'], 'w');
soltype = 4;
for i = dims
    N = 2^i;
    M = N-50;
    F = randn(M,N);
    % F = triu(F,3) + tril(F,3);
    if soltype == 1
        f = zeros(N,1);
        f(N/4) = 1;
        M = N-200;
        b = fft(f);
        b = b(1:M);
    elseif soltype == 2        
        f = ones(N,1);
        f(N/2+1:N) = 2;
        M = N-20;
        b = fft(f);
        b = b(1:M);
    elseif soltype == 3
        f = zeros(N,1);
        f(N/4) = 1;
        b = F*f;
        % b = b(1:M);
    elseif soltype == 4
        f = ones(N,1);
        f(N/2+1:N) = 2;
        % f(N/4+1:N/2) = 1.001;
        b = F*f;
    end
    % hf = @(x)x(1:M);
    hfM = @(x)x(1:M);
    % ft = [ones(M,1);zeros(N/2,1)];
    ftM = [ones(M,1);zeros(N-M,1)];
    f0s = [f 2*f randn(N,1)];
    for mu = mus
        lmus = [mu/8 mu/4 mu/2 mu mu*2 mu*4 mu*8]
        for lambda = lmus
            data = [];
            for kk = 1:size(f0s,2)
                f0 = f0s(:,kk);
                Inv = [ones(M,1)/(mu*N+lambda);ones(N-M,1)/lambda];
                InvF = [ones(M,1)/(mu*N-2*lambda);ones(N-M,1)/(-2*lambda)];
                if soltype == 3
                    [L,~] = chol(mu*F'*F+lambda*eye(N));
                elseif soltype == 2
                    BB = mu*N*ifft(ftM.*fft(eye(N)))+lambda*laplacianp(eye(N), 1/N^2);
                    img = imag(diag(BB));
                    [L,D] = ldl(BB-1i*diag(img));
                    % [Fac,HODLR] = HODLR_construction_hmt( N, @(x)mu*N*ifft(fft(x).*ftM)+lambda*laplacianp(x, 1/N^2), 1E-12, NoFIle, 32, 64,64);
                    % [G] = RSS_ldl(Fac,1E-12,NoFIle);
                elseif soltype == 4
                    [L,D] = ldl(mu*F'*F+lambda*laplacianp(eye(N),1/N^2));
                    % [Fac,HODLR] = HODLR_construction_hmt( N, @(x)mu*F'*F*x+lambda*laplacianp(x,1/N^2), 1E-12, NoFIle, 32, 64,64);
                    % [G] = RSS_ldl(Fac,1E-12,NoFIle);
                end
                % [x, iter] = SplitBregman('L1', mu, lambda, 1/N, @(x)N*ifft(x), @(x)mu*N*ifft(fft(x))+lambda*x, b, @(x)x, @(x)x, maxit1, maxit2, tol, @(x)ifft(fft(x))/(N*mu+lambda), f0, f0, zeros(N,1));
                % [x, iter] = SplitBregman('L1', mu, lambda, 1/N, @(x)N*ifft([x;zeros(N-M,1)]), @(x)mu*N*ifft(fft(x).*ftM)+lambda*x, b, @(x)x, @(x)x, maxit1, maxit2, tol, @(x)ifft(fft(x).*Inv), f0, f0, zeros(N,1));
                % [x, iter] = SplitBregman('L1', mu, lambda, 1/N, @(x)F'*x, @(x)mu*F'*F*x+lambda*x, b, @(x)x, @(x)x, maxit1, maxit2, tol, @(x)L\(L'\x), f0, f0, zeros(N,1));
                % [x, iter] = SplitBregman('TV-L1', mu, lambda, 1/N, @(x)N*ifft(x), @(x)mu*N*ifft(fft(x))+lambda*laplacianp(x, 1/N^2), b, @(x)gradientxp(x,1/N), @(x)x, maxit1, maxit2, tol, @(x)ifft(fft(x))/mu/N, f0, f0,zeros(N,1));
                % [x, iter] = SplitBregman('TV-L1', mu, lambda, 1/N, @(x)N*ifft([x;zeros(N-M,1)]), @(x)mu*N*ifft(fft(x).*ftM)+lambda*laplacianp(x, 1/N^2), b, @(x)gradientxp(x,1/N), @(x)gradientxp(x,1/N), maxit1, maxit2, tol, @(x)L'\(D\(L\x)), f0, f0,zeros(N,1));
                [x, iter] = SplitBregman('TV-L1', mu, lambda, 1/N, @(x)F'*x, @(x)mu*F'*F*x+lambda*laplacianp(x,1/N^2), b, @(x)gradientxp(x,1/N), @(x)gradientxp(x,1/N), maxit1, maxit2, tol, @(x)L'\(D\(L\x)), f0, f0,zeros(N,1));
                
                % [x, iter] = SplitBregman('TV-L1', mu, lambda, 1/N, @(x)N*ifft(x), @(x)mu*N*ifft(fft(x))+lambda*laplacianp(x, 1/N^2).*gradientxp(x,1/N), b, @(x)laplacianp(x,1/N^2), maxit1, maxit2, tol, @(x)ifft(fft(x))/mu/N, f, f,zeros(N,1));
                % [x, iter] = SplitBregman('TV-L1', mu, lambda, 1/N, @(x)N*ifft(x), @(x)mu*N*ifft(fft(x))+lambda*sum(laplacianp(x, 1/N^2)), b, @(x)gradientxpT(x,1/N), maxit1, maxit2, tol, @(x)ifft(fft(x))/mu/N, f, f,zeros(N,1));
                % [x, iter] = SplitBregman('TV-L1', mu, lambda, 1/N, @(x)N*ifft(x), @(x)mu*N*ifft(fft(x))+lambda*gradientxpT(gradientxp(x,1/N),1/N), b, @(x)gradientxpT(x,1/N), maxit1, maxit2, tol, @(x)ifft(fft(x))/mu/N, f, f,zeros(N,1));
                if soltype < 3
                    y = fft(x);
                    e1 = norm(b-y(1:M))/norm(b);
                else
                    y = F*x;
                    e1 = norm(b-y)/norm(b);
                end
                gx = gradientxp(x,1/N);
                gx(N/4-5:N/4+5)'
                e2 = norm(x-f)/norm(f);
                % e = norm(f-x)/norm(f);
                data = [data iter(1) iter(2) e1 e2];
            end
            fprintf(OutPutFile,'\\hline \n')
            fprintf(OutPutFile,'($2^{%d }$, $2^{%d }$) & %d & %d & %5.1e & %5.1e & %d & %d & %5.1e & %5.1e & %d & %d & %5.1e & %5.1e \\\\ \n', log2(mu), log2(lambda), data(1), data(2), data(3), data(4), data(5), data(6), data(7), data(8), data(9), data(10), data(11), data(12))
        end
    end
    % fprintf('results %10.4e / %d / %d \n', e, iter(1), iter(2))

end