function v = laplacianp(f, h)
    [n,m] = size(f);
    v = zeros(size(f));
    v(1,:) = (f(n,:) - 2*f(1,:) + f(2,:))/h;
    v(n,:) = (f(n-1,:) - 2*f(n,:) + f(1,:))/h;
    for i = 2:n-1
        v(i,:) = (f(i-1,:) - 2*f(i,:) + f(i+1,:))/h;
    end
end