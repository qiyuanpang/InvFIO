function a0 = linesearch(fun, x0, p, m, tau, c)
    a0 = 5;
    t = -c*m;
    while fun(x0) - fun(x0+a0*p) < a0*t
        a0 = a0*tau;
    end
end