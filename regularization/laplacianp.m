function v = laplacianp(f, h)
    [n,m] = size(f);
    v = zeros(size(f));
    v(1,:) = (f(n,:) - 2*f(1,:) + f(2,:))/h/h;
    v(n,:) = (f(n-1,:) - 2*f(n,:) + f(1,:))/h/h;
    v(2:n-1,:) = (f(1:n-2,:) - 2*f(2:n-1,:) + f(3:n,:))/h/h;
end
