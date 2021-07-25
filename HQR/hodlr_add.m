function H = hodlr_add(H1, H2, tol_or_rank)
    if H1.lvl == 0
        H = struct('A11', H1.A11+H2.A11, 'A12', [], 'A21', [], 'A22', [], 'lvl', H1.lvl, 'm', H1.m, 'n', H1.n, 'maxrk', 0);
    else
        H = struct('A11', [], 'A12', [], 'A21', [], 'A22', [], 'lvl', H1.lvl, 'm', H1.m, 'n', H1.n, 'maxrk', 0);
        H.A12 = struct('Q', [], 'R', [], 'm', H1.A12.m, 'n', H1.A12.n, 'r', 0);
        [H.A12.Q, H.A12.R] = recompression([H1.A12.Q H2.A12.Q], [H1.A12.R; H2.A12.R], tol_or_rank);
        H.A12.r = size(H.A12.Q, 2);
        H.A21 = struct('Q', [], 'R', [], 'm', H1.A21.m, 'n', H1.A21.n, 'r', 0);
        [H.A21.Q, H.A21.R] = recompression([H1.A21.Q H2.A21.Q], [H1.A21.R; H2.A21.R], tol_or_rank);
        H.A21.r = size(H.A21.Q, 2);
        H.maxrk = max(H.A12.r, H.A21.r);
        H.A11 = hodlr_add(H1.A11, H2.A11, tol_or_rank);
        H.A22 = hodlr_add(H1.A22, H2.A22, tol_or_rank);
    end
end