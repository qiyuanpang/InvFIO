function x = hodlr_tri2_sol(F1, F2, b)
    if F1.lvl == 0
        x = (F1.A11+F2.A11)\b;
    else
        x = zeros(size(b));
        n = size(b, 1);
        n2 = floor(n/2);
        if isempty(F1.A21)
            x(n2+1:n,:) = hodlr_tri2_sol(F1.A22, F2.A22, b(n2+1:n,:));
            c = F1.A12.Q*(F1.A12.R*x(n2+1:n,:))+F2.A12.Q*(F2.A12.R*x(n2+1:n,:));
            x(1:n2,:) = hodlr_tri2_sol(F1.A11, F2.A11, b(1:n2,:)-c);
        else if isempty(F1.A12)
            x(1:n2,:) = hodlr_tri2_sol(F1.A11, F2.A11, b(1:n2));
            c = F1.A21.Q*(F1.A21.R*x(1:n2,:))+F2.A21.Q*(F2.A21.R*x(1:n2,:));
            x(n2+1:n,:) = hodlr_tri2_sol(F1.A22, F2.A22, b(n2+1:n,:)-c);
        end
    end
end
