function y = hodlr_adj_apply(F, v)
    if F.lvl > 0
        n = size(v ,1);
        y = zeros(size(v));
        n2 = floor(n/2);
        if isempty(F.A21)
            y(1:n2, :) = hodlr_adj_apply(F.A11, v(1:n2, :));
        else
            y(1:n2, :) = hodlr_adj_apply(F.A11, v(1:n2, :)) + F.A21.R'*(F.A21.Q'*v(n2+1:n, :));
        end
        if isempty(F.A12)
            y(n2+1:n, :) = hodlr_adj_apply(F.A22, v(n2+1:n, :));
        else
            y(n2+1:n, :) = hodlr_adj_apply(F.A22, v(n2+1:n, :)) + F.A12.R'*(F.A12.Q'*v(1:n2, :));
        end
    else
        y = F.A11'*v;
    end
end