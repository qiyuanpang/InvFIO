function v = laplacian(f, h)
    n = size(f,1);
    v = zeros(size(f));
    v(1,:) = (f(3,:)+f(1,:)-2*f(2,:))/h/h;
    v(n,:) = (f(n,:)+f(n-2,:)-2*f(n-1,:))/h/h;
    v(2:n-1,:) = (f(3:n,:)+f(1:n-2,:)-2*f(2:n-1,:))/h/h;
end