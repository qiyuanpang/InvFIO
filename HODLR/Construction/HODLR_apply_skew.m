function [ res ] = Hmatrix_1d_apply_skew( F,x )
    %HMATRIX_1D_APPLY
    % Function to apply H_matrix representation F to a vector x.
    res = zeros(size(x));
    for nbk = 1:F(1).nb
       res(F(nbk).x,:)= res(F(nbk).x,:) + F(nbk).U*(F(nbk).S*(F(nbk).V'*x(F(nbk).y,:)));
       res(F(nbk).y,:)= res(F(nbk).y,:) - F(nbk).V*(F(nbk).S'*(F(nbk).U'*x(F(nbk).x,:)));
    end
    
    nd = F(1).nd;
    if ~isempty(nd)
        for nbk = nbk+1:F(1).nd
            res(F(nbk).x,:)= res(F(nbk).x,:) + F(nbk).U*x(F(nbk).x,:);
        end
    end
    
end