function [L, U] = HODLR_ldu(H, lambda)
    if H.lvl > 0
        L = struct('A11', [], 'A12', [], 'A21', [], 'A22', [], 'lvl', 0, 'm', 0, 'n', 0, 'maxrk', 0);
        U = struct('A11', [], 'A12', [], 'A21', [], 'A22', [], 'lvl', 0, 'm', 0, 'n', 0, 'maxrk', 0);
        L.A21 = H.A21;
        U.A12 = H.A12;
        L.m = H.m;
        U.m = H.m;
        L.n = H.n;
        U.n = H.n;
        L.maxrk = H.maxrk;
        U.maxrk = H.maxrk;
        L.lvl = H.lvl;
        U.lvl = H.lvl;
        [L.A11, U.A11] = HODLR_ldu(H.A11, lambda);
        [L.A22, U.A22] = HODLR_ldu(H.A22, lambda);
    else
        L = struct('A11', tril(H.A11)+eye(size(H.A11))*lambda, 'A12', [], 'A21', [], 'A22', [], 'lvl', 0, 'm', 0, 'n', 0, 'maxrk', 0);
        U = struct('A11', triu(H.A11,1), 'A12', [], 'A21', [], 'A22', [], 'lvl', 0, 'm', 0, 'n', 0, 'maxrk', 0);
        L.m = H.m;
        U.m = H.m;
        L.n = H.n;
        U.n = H.n;
        L.maxrk = H.maxrk;
        U.maxrk = H.maxrk;
        L.lvl = H.lvl;
        U.lvl = H.lvl;
    end
end