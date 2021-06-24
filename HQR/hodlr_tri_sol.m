function x = hodlr_tri_sol(F, b)
    if F.lvl == 0
        x = F.A11\b;
    else
        x = zeros(size(b));
        n = size(b, 1);
        n2 = floor(n/2);
        x(n2+1:n,:) = hodlr_tri_sol(F.A22, b(n2+1:n,:));
        c = F.A12.Q*(F.A12.R*x(n2+1:n,:));
        x(1:n2,:) = hodlr_tri_sol(F.A11, b(1:n2,:)-c);
    end
end