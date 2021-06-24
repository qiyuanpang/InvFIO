function [Y T R] = householderqr(A)
    [m n] = size(A);
    Y = zeros(m, n);
    T = eye(n)*(-2);
    R = A;
    if m < n
        fprintf("Input matrix should be a tall skinny or square matrix!")
        return;
    end
    for k = 1:n
        x = R(k:m, k);
        x(1) = sign(x(1))*norm(x) + x(1);
        x = x/norm(x);
        Y(k:m,k) = x;
        R(k:m,k:n) = R(k:m,k:n) - 2*x*(x'*R(k:m,k:n));
        if k > 1
            T(1:k-1,k) = -2*T(1:k-1,1:k-1)*(Y(k:m,1:k-1)'*x);
        end
    end
    R = R(1:n,1:n);
    T = -T;
end
