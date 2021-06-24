function [U,D,V] = RSVD( A, tol, r )
%RSVD Summary of this function goes here
%   Detailed explanation goes here
[U,D,V] = svd(A,'econ');

if D(1,1)<tol
    U=[];
    V=[];
    D=[];
else
    idx = max(find(find(diag(D)>tol*D(1,1))<=r));
    U = U(:,1:idx);
    V = V(:,1:idx);
    D = D(1:idx,1:idx);
end

end

