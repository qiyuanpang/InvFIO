function [u, iters] = SplitBregman_TV2D(mu, lambda, funAT, funH, f, opR0, opR1, Phi0, Phi1, maxit1, maxit2, tol, funM, u0, d0, b0, d1,b1)
    N = size(f,1);
    % [L, U] = HODLR_ldu(H, 0);
    % [L1, U1] = HODLR_laplacian(H, -lambda/mu/h/h);
    % HODLR_laplacian_modA11(L1, U1, -lambda/mu/h/h);
    % HODLR_laplacian_modA22(L1, U1, -lambda/mu/h/h);
    % fprintf('cond %10.4e %10.4e \n', cond(funH(eye(N))), cond(funH(funM(eye(N)))))
    

    iter1 = 0;
    for iter = 1:maxit1
        %[u, iter1] = GaussSeidel(@(x)1/mu*hodlr_tri2_sol(L,L1,x), @(x)mu*hodlr_apply(U,x)+mu*hodlr_apply(U1,x), mu*funAT(f)+lambda*opR(d0-b0), u0, maxit2, tol);
        for k = 1:1
            [u, flag, relres, iter2] = pcg(funH, mu*funAT(f)+lambda*opR0(d0-b0)+lambda*opR1(d1-b1), 1E-14, maxit2, funM, @(x)x, u0);
            % norm(funH(u0)-funAT(f))/norm(funAT(f))
            % norm(funM(funH(u0))-u0)/norm(u0)
            % [u, flag, relres, iter2] = gmres(funH, mu*funAT(f)+lambda*opR(d0-b0), 5, 1E-13, maxit2, funM, @(x)x, u0);
            % norm(funH(u) - mu*funAT(f)-lambda*opR(d0-b0))/norm(mu*funAT(f)+lambda*opR(d0-b0))
            % cond(funM(funH(eye(N))))
            % max(max(abs(u)))
            % fprintf('pcg %10.4e / %d / %10.4e \n', norm(u-u0), flag, relres)
            iter1 = max(iter1, iter2);
            if norm(u-u0) < tol
                break
            end
            u0 = u;
            d0 = shrink(Phi0(u0)+b0, 1/lambda);
            d1 = shrink(Phi1(u0)+b1, 1/lambda);
        end
        b0 = b0 + Phi0(u0) - d0;
        b1 = b1 + Phi1(u0) - d1;
    end
    iters = [iter, iter1];
    
end
