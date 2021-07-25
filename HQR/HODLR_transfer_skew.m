function F = HODLR_transfer_skew(H, maxlvl, lvl, root)
    if lvl > 0
        d = maxlvl - lvl;
        prev1 = 2^d - 1;
        prev2 = root - prev1 -1;
        prev = 2^(d+1) - 1 + 2*prev2;
        F = struct('A11', [], 'A12', [], 'A21', [], 'A22', [], 'lvl', lvl, 'm', 0, 'n', 0, 'maxrk', 0);
        F.A21 = struct('Q', H(root).U, 'R', H(root).S*H(root).V',  'lvl', lvl, 'm', size(H(root).U, 1), 'n', size(H(root).V, 1), 'r', size(H(root).U, 2));
        F.A12 = struct('Q', H(root).V, 'R', -H(root).S'*H(root).U',  'lvl', lvl, 'm', size(H(root).V, 1), 'n', size(H(root).U, 1), 'r', size(H(root).V, 2));
        F.m = F.A21.m + F.A12.m;
        F.n = F.A21.n + F.A12.n;
        F.A11 = HODLR_transfer_skew(H, maxlvl, lvl-1, prev+1);
        F.A22 = HODLR_transfer_skew(H, maxlvl, lvl-1, prev+2);
        F.maxrk = max([F.A12.r F.A21.r F.A11.maxrk F.A22.maxrk]);
    else
        F = struct('A11', H(root).U, 'A12', [], 'A21', [], 'A22', [], 'lvl', lvl, 'm', size(H(root).U, 1), 'n', size(H(root).U,2), 'maxrk', 0);
    end
end