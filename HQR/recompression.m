function [Q, R] = recompression(Al, Ar, tol)
    if size(Al,2) == 1
        Q = Al;
        R = Ar;
    else
        [Q1, R1] = qr(Al, 0);
        [Q2, R2] = qr(Ar', 0);
        [U, S, V] = svd(R1*R2', 0);
        diagS = diag(S);
        k = nnz(abs(diagS) > tol);
        Q = Q1*U(:, 1:k);
        R = diagS(1:k).*(V(:,1:k)'*Q2');
    end
end
