function [x] = RSS_apply_lu(F,x)
%RSS_APPLY summary of this function goes here
%   Detailed explanation goes here

for n_f=1:F.nblocks-1
    rd = F.factors(n_f).rd;
    sk = F.factors(n_f).sk;
    if ~isempty(rd)
       x(sk,:) = x(sk,:)+F.factors(n_f).T*x(rd,:); 
       x(rd,:) = F.factors(n_f).U*x(rd,:);
       x(rd,:) = x(rd,:) + F.factors(n_f).C2*x(sk,:);
    end
end

%Apply last matrix level
sk = F.factors(F.nblocks).sk;
x(sk,:) = F.factors(F.nblocks).L*F.factors(F.nblocks).U*x(sk,:);

for n_f=F.nblocks-1:-1:1
    rd = F.factors(n_f).rd;
    sk = F.factors(n_f).sk;
    if ~isempty(rd)
       x(sk,:) = x(sk,:)+F.factors(n_f).C1*x(rd,:);
       x(rd,:) = F.factors(n_f).L*x(rd,:);
       x(rd,:) = x(rd,:) + F.factors(n_f).T'*x(sk,:);
    end
end


end

