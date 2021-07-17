function HODLR_laplacian_modA11(L, U, h)
    if L.lvl > 0
        HODLR_laplacian_modA11(L.A11, U.A11, h);
    else
        L.A11(1) = h;
        U.A11(1,2) = -2*h;
        U.A11(1,3) = h;
    end
end