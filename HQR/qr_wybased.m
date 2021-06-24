function [Y T R] = qr_wybased(A)
    [m n] = size(A);
    if m < n
        fprintf("Input matrix should be a tall skinny or square matrix!")
        return;
    end
    %% if n == 1, return (I - 2 v v^T)*(rho 0)^{T}
    if n < max(8, floor(log2(n+1)))
        [Y T R] = householderqr(A);
        %x = A(:);
        %v = sign(x(1))*norm(x)*[1; zeros(m-1, 1)] + x;
        %Y = v / norm(v);
        %T = [2];
        %R = [x(1)-2*Y(1)*Y'*x];
    else
        n1 = floor(n/2);
        [Y1 T1 R1] = qr_wybased(A(:, 1:n1));
        A(:, n1+1:n) = A(:, n1+1:n) - Y1*(T1'*(Y1'*A(:, n1+1:n)));
        [Y2 T2 R2] = qr_wybased(A(n1+1:m, n1+1:n));
        %Y2 = [zeros(n1, n-n1); Y2];
        T12 = -(T1*Y1(n1+1:m, :)')*Y2*T2;
        Y2 = [zeros(n1, n-n1); Y2];
        Y = [Y1 Y2];
        T = [T1 T12; zeros(n-n1, n1) T2];
        R = [R1 A(1:n1, n1+1:n); zeros(n-n1, n1) R2];
    end

end
