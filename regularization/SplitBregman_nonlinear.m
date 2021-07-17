function [u0, iter] = SplitBregman_nonlinear(funA, funAT, f, Phi, Dphi, lambda, tol, maxit, u0, d0, b0)
    switch nargin
        case 7
            maxit = 50;
            u0 = zeros(size(f));
            d0 = zeros(size(f));
            b0 = zeros(size(f));
        case 8
            u0 = zeros(size(f));
            d0 = zeros(size(f));
            b0 = zeros(size(f));
        case 9
            d0 = zeros(size(f));
            b0 = zeros(size(f));
        case 10
            b0 = zeros(size(f));
    end
    for iter = 1:maxit
        fun = @(x)norm(funA(x)-f)^2/norm(f)^2 + lambda/2*norm(d0 - Phi(x) - b0)^2;
        Dfun = @(x)(2*funAT(funA(x))-2*funAT(f))/norm(f)^2 + lambda*Dphi(x)*(Phi(x)-d0+b0);
        % dx = zeros(size(u0));
        % dx(1) = 1E-8;
        % Df = Dfun(u0);
        % norm((fun(u0+dx)-fun(u0))/norm(dx) - Df(1))
        u1 = u0;
        if iter == 1
            fprintf('initial loss: %10.4e \n', fun(u0))
        end
        [u0, iter1] = NonlinearCG(fun, Dfun, tol, u1, maxit);
        fprintf('iter: %d / %10.4e \n', iter, fun(u0))
        % norm(d0 - shrink(Phi(u0)+b0, 1/lambda))
        d0 = shrink(Phi(u0)+b0, 1/lambda);
        b0 = b0 + Phi(u0) - d0;
        % norm(u0 - u1)
        if norm(u0-u1) < tol
            break
        end
    end

end