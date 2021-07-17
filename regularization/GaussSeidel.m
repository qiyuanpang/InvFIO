function [x, iter] = GaussSeidel(funInvL, funU, b, x0, maxit, tol)
    for iter = 1:maxit
        x = funInvL(b) - funInvL(funU(x0));
        if norm(x-x0) < tol
            break;
        x0 = x;
    end
end
