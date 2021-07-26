function v = gradientxpT(f, h)
    n = size(f,1);
    v = zeros(size(f));
    v(1,:) = -(f(2,:)-f(n,:))/2/h;
    v(n,:) = -(f(1,:)-f(n-1,:))/2/h;
    v(2:n-1) = -(f(3:n)-f(1:n-2))/2/h;
end