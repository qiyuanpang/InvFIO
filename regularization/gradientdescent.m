function [x1, iter] = gradientdescent(Dfun, x0, tol)
    iter = 1;
    x1 = x0 - 0.1*Dfun(x0);
    while norm(x1-x0) >= tol
        Df1 = Dfun(x1);
        Df0 = Dfun(x0);
        gamman = abs((x1-x0)'*(Df1-Df0))/norm(Df1-Df0)^2;
        x0 = x1;
        x1 = x1 - gamman*Df1;
        iter = iter + 1;
    end
end