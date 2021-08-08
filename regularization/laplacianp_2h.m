function v = laplacianp_2h(f, h)
    [n,m] = size(f);
    v = zeros(size(f));
    v(1,:) = (f(n-1,:) - 2*f(1,:) + f(3,:))/4/h/h;
    v(2,:) = (f(n,:) - 2*f(2,:) + f(4,:))/4/h/h;
    v(n,:) = (f(n-2,:) - 2*f(n,:) + f(2,:))/4/h/h;
    v(n-1,:) = (f(n-3,:) - 2*f(n-1,:) + f(1,:))/4/h/h;
    v(3:n-2,:) = (f(1:n-4,:) - 2*f(3:n-2,:) + f(5:n,:))/4/h/h;
end
