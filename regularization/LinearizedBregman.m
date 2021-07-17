function [u0, iter] = LinearizedBregman(funA, funAT, f, delta, mu, tol, maxit, u0, v0)
    switch nargin
        case 6
            maxit = 50;
            u0 = zeros(size(f));
            v0 = zeros(size(f));
        case 7
            u0 = zeros(size(f));
            v0 = zeros(size(f));
        case 8
            v0 = zeros(size(f));
        otherwise
            error("Input parameters are too many or too few!")
    end
    for iter = 1:maxit
        err = norm(f-funA(u0))/norm(f);
        if err < tol
            break
        end
        v0 = v0 + funAT(f - funA(u0));
        u0 = delta*shrink(v0, 1/mu);
    end
        
end