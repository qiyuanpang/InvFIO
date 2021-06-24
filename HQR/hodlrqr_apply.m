function y = hodlrqr_apply(Y, T, R, x)
    y = hodlr_apply(R, x);
    z = hodlr_adj_apply(Y,y);
    y = y - hodlr_apply(Y, hodlr_apply(T, z));
end