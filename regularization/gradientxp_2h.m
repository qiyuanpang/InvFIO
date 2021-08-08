function v = gradientxp_2h(f, h)
    n = size(f,1);
    v = zeros(size(f));
    v(1,:) = (f(3,:)-f(n-1,:))/4/h;
    v(2,:) = (f(4,:)-f(n,:))/4/h;
    v(n,:) = (f(2,:)-f(n-2,:))/4/h;
    v(n-1,:) = (f(1,:)-f(n-3,:))/4/h;
    v(3:n-2) = (f(5:n,:)-f(1:n-4,:))/4/h;
end
