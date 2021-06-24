function [ res ] = HODLR_RSS_apply( F,x )
%HMATRIX_1D_APPLY when factors are stored in RSS way.
% Function to apply HODLR representation F to a vector x.
nb = 2^(F.nlvl-1); 
res = zeros(size(x));
for lvl = 1:F.nlvl-1
    for nbk = 1:2:2^lvl
       blk_size = N/2^lvl;
       res(F(nbk).x,:)= res(F(nbk).x,:) + F(nbk).U*F(nbk).V'*x(F(nbk).y,:);
       res(F(nbk).y,:)= res(F(nbk).y,:) + F(nbk).V*F(nbk).U'*x(F(nbk).x,:);
    end
end
    for nbk_n = 1:nb 
      res(F.factors(nbk).slf,:)= res(F.factors(nbk).slf,:) + F.factors(nbk).V*F.factors(nbk_n).V'*x(F.factors(nbk_n).slf,:);
    end
    
for nbk = 1:nb
   res(F.factors(nbk).slf,:)= res(F.factors(nbk).slf,:) + F.factors(nbk).A_blkdiag*x(F.factors(nbk).slf,:);
end

end


