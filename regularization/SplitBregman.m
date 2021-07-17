function [u, iter] = SplitBregman(regm, mu, lambda, h, funAT, H, f, opR, maxit1, maxit2, tol, funM, u0, d0, b0)
    if strcmp(regm, 'L1')
        % [L, U] = HODLR_ldu(H, lambda/mu);
        for iter = 1:maxit1
            % [u, iter1] = GaussSeidel(@(x)1/mu*hodlr_tri_sol(L,x), @(x)mu*hodlr_apply(U,x), mu*funAT(f)+lambda*opR(d0-b0), u0, maxit, tol);
            [u, flag, relres, iter1] = pcg(@(x)mu*hodlr_apply(H,x)+lambda*x, mu*funAT(f)+lambda*(d0-b0), tol, maxit2, @(x)1/mu*funM(x), @(x)x, u0);
            if norm(u-u0) < tol
                break
            end
            u0 = u;
            d0 = shrink(opR(u0)+b0, 1/lambda);
            b0 = b0 + opR(u0) - d0;
        end
    else if strcmp(regm, 'TV-L1')
        N = size(f,1);
        [L, U] = HODLR_ldu(H, 0);
        [L1, U1] = HODLR_laplacian(H, -lambda/mu/h/h);
        HODLR_laplacian_modA11(L1, U1, -lambda/mu/h/h);
        HODLR_laplacian_modA22(L1, U1, -lambda/mu/h/h);
        for iter = 1:maxit1
            % [u, iter1] = GaussSeidel(@(x)1/mu*hodlr_tri2_sol(L,L1,x), @(x)mu*hodlr_apply(U,x)+mu*hodlr_apply(U1,x), mu*funAT(f)+lambda*opR(d0-b0), u0, maxit2, tol);
            [u, flag, relres, iter1] = pcg(@(x)mu*(hodlr_apply(H,x)+hodlr_apply(L1,x)+hodlr_apply(U1,x)), mu*funAT(f)+lambda*opR(d0-b0), tol, maxit2, @(x)funM(x), @(x)x, u0);
            if norm(u-u0) < tol
                break
            end
            u0 = u;
            d0 = shrink(opR(u0)+b0, 1/lambda);
            b0 = b0 + opR(u0) - d0;
        end
    end
end