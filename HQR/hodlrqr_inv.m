function x = hodlrqr_inv(Y, T, R, b)
    z = b - hodlr_apply(Y, hodlr_adj_apply(T, hodlr_adj_apply(Y, b)));
    x = hodlr_tri_sol(R, z);
end
