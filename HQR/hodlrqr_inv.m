function x = hodlrqr_inv(Y, T, R, b)
    fun = @(b)(b - hodlr_apply(Y, hodlr_adj_apply(T, hodlr_adj_apply(Y, b))));
    z = fun(b);
    %z = gmres(@(x)(x - hodlr_apply(Y, hodlr_apply(T, hodlr_adj_apply(Y, x)))), b, 3, 1E-10, 6, fun, @(x)(x), fun(b));
    x = hodlr_tri_sol(R, z);
end
