function F = hodlr_sub_lr(H, Q, R, tol_or_rank)
    if H.lvl == 0
        A11 = H.A11 - Q*R;
        F = struct('A11', A11, 'A12', [], 'A21', [], 'A22', [], 'lvl', H.lvl, 'm', H.m, 'n', H.n, 'maxrk', 0);
    else
        F = struct('A11', [], 'A12', [], 'A21', [], 'A22', [], 'lvl', H.lvl, 'm', H.m, 'n', H.n, 'maxrk', 0);
        m1 = H.A11.m;
        n1 = H.A11.n;
        m2 = H.A22.m;
        n2 = H.A22.n;
        F.A12 = struct('Q', [], 'R', [], 'm', m1, 'n', n2, 'r', 0);
        [F.A12.Q, F.A12.R] = recompression([H.A12.Q -Q(1:m1,:)], [H.A12.R; R(:,n1+1:n1+n2)], tol_or_rank);
        F.A12.r = size(F.A12.Q, 2);
        F.A21 = struct('Q', [], 'R', [], 'm', m2, 'n', n1, 'r', 0);
        [F.A21.Q, F.A21.R] = recompression([H.A21.Q -Q(m1+1:m1+m2,:)], [H.A21.R; R(:,1:n1)], tol_or_rank);
        F.A21.r = size(F.A21.Q, 2);
        F.A11 = hodlr_sub_lr(H.A11, Q(1:m1,:), R(:,1:n1), tol_or_rank);
        F.A22 = hodlr_sub_lr(H.A22, Q(m1+1:m1+m2,:), R(:,n1+1:n1+n2), tol_or_rank);
        F.maxrk = max([F.A12.r F.A21.r F.A11.maxrk F.A22.maxrk]);
    end
end