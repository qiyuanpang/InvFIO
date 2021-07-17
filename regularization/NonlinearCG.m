function [x, iter] = NonlinearCG(fun, Dfun, tol, x0, maxit)
    switch nargin 
        case 4
            maxit = 20;
    end
    dx0 = -Dfun(x0);
    % alpha0 = linesearch(fun, x0, dx0, dx0'*dx0, 0.5, 0.5);
    [alpha0, iter0] = gradientdescent(@(t)real(Dfun(x0+t*dx0)'*dx0), 10.0, 1E-10);
    s = dx0;
    x = x0 + alpha0*dx0;
    for iter = 1:maxit
        if norm(x-x0) < tol
            break
        end
        dx = -Dfun(x);
        betaPR = (dx'*(dx - dx0))/(dx0'*dx0);
        % betan = max(0, abs(betaPR)*sign(real(betaPR)));
        betan = max(0, abs(betaPR));
        s = dx + betan*s;
        dx0 = dx;
        % alpha0 = linesearch(fun, x, s, dx0'*s, 0.5, 0.5);
        [alpha0, iter0] = gradientdescent(@(t)real(Dfun(x+t*s)'*s), 10.0, 1E-10);
        % fprintf('Line: %10.4e \n', fun(x)-fun(x+alpha0*s))
        x0 = x;
        x = x + alpha0*s;
    end
end