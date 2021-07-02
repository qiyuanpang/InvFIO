function [ x ] = RSS_inv_lu(F,x)
%RSS_INV Summary of this function goes here
%   Detailed explanation goes here

for n_f=1:F.nblocks-1
    rd = F.factors(n_f).rd;
    sk = F.factors(n_f).sk;
    if ~isempty(rd)
       x(rd,:) = x(rd,:)-F.factors(n_f).T'*x(sk,:); 
       x(rd,:) = F.factors(n_f).L\x(rd,:);
       x(sk,:) = x(sk,:) - F.factors(n_f).C1*x(rd,:);
    end
end

%Invert last level matrix
sk = F.factors(F.nblocks).sk;
x(sk,:) = F.factors(F.nblocks).U\(F.factors(F.nblocks).L\x(sk,:));

for n_f=F.nblocks:-1:1
    rd = F.factors(n_f).rd;
    sk = F.factors(n_f).sk;
    if ~isempty(rd)
       x(rd,:) = x(rd,:)-F.factors(n_f).C2*x(sk,:);
       x(rd,:) = F.factors(n_f).U\x(rd,:);
       x(sk,:) = x(sk,:) - F.factors(n_f).T*x(rd,:);
    end
end

end

