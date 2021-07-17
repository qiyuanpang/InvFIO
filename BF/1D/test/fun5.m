function res = fun5(N,x,k)

    x = (x-1)/N;
    k = k-1-N/2;
    xk = x*k';
    sx = (2 + sin(2*pi*x))/5.94;
    tmp = (2*pi)* (xk + sx*abs(k'));
    res = complex(cos(tmp),sin(tmp));
    
end
