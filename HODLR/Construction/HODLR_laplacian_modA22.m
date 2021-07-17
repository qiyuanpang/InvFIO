function HODLR_laplacian_modA22(L, U, h)
    if L.lvl > 0
        HODLR_laplacian_modA22(L.A22, U.A22, h);
    else
        L.A11(L.m,L.n) = h;
        L.A11(L.m,L.n-1) = -2*h;
        L.A11(L.m,L.n-2) = h;
    end
end