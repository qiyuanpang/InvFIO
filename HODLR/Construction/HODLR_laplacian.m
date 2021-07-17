function [L, U] = HODLR_laplacian(H, h)
    if H.lvl > 0
        L = struct('A11', [], 'A12', [], 'A21', [], 'A22', [], 'lvl', H.lvl, 'm', H.m, 'n', H.n, 'maxrk', H.maxrk);
        U = struct('A11', [], 'A12', [], 'A21', [], 'A22', [], 'lvl', H.lvl, 'm', H.m, 'n', H.n, 'maxrk', H.maxrk);
        L.A21 = struct('Q', zeros(H.A21.m,1), 'R', zeros(1,H.A21.n),  'lvl', H.A21.lvl, 'm', H.A21.m, 'n', H.A21.n, 'r', H.A21.r);
        L.A21.Q(1) = 1;
        L.A21.R(end) = h;
        U.A12 = struct('Q', zeros(H.A12.m,1), 'R', zeros(1,H.A12.n),  'lvl', H.A12.lvl, 'm', H.A12.m, 'n', H.A12.n, 'r', H.A12.r);
        U.A12.Q(end) = 1;
        U.A12.R(1) = h;
        [L.A11, U.A11] = HODLR_laplacian(H.A11, h);
        [L.A22, U.A22] = HODLR_laplacian(H.A22, h);
    else
        L = struct('A11', [], 'A12', [], 'A21', [], 'A22', [], 'lvl', H.lvl, 'm', H.m, 'n', H.n, 'maxrk', H.maxrk);
        L.A11 = eye(L.m,L.n)*(-2)*h;
        L.A11(2:L.m,1:L.n-1) = h;
        U = struct('A11', [], 'A12', [], 'A21', [], 'A22', [], 'lvl', H.lvl, 'm', H.m, 'n', H.n, 'maxrk', H.maxrk);
        U.A11 = zeros(U.m, U.n);
        U.A11(1:U.m-1,2:U.n) = h;
    end
end