function y = hodlr_apply(F, v)
    if F.lvl > 0
        n = size(v ,1);
        y = zeros(size(v));
        n2 = floor(n/2);
        if isempty(F.A12)
            y(1:n2, :) = hodlr_apply(F.A11, v(1:n2, :));
        else
            y(1:n2, :) = hodlr_apply(F.A11, v(1:n2, :)) + F.A12.Q*(F.A12.R*v(n2+1:n, :));
        end
        if isempty(F.A21)
            y(n2+1:n, :) = hodlr_apply(F.A22, v(n2+1:n, :));
        else
            y(n2+1:n, :) = hodlr_apply(F.A22, v(n2+1:n, :)) + F.A21.Q*(F.A21.R*v(1:n2, :));
        end
    else
        y = F.A11*v;
    end
end