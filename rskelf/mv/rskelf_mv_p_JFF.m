% RSKELF_MV_P  Dispatch for RSKELF_MV with F.SYMM = 'P'.
%
%    See also RSKELF, RSKELF_MV.

function Y = rskelf_mv_p_JFF(F,X)

  % initialize
  n = F.lvp(end);
  Y = X;

  % upward sweep
  for i = 1:n
    sk = F.factors(i).sk;
    rd = F.factors(i).rd;
    slf = F.factors(i).slf; %JFF addition
    Y(slf,:) = F.factors(i).D*Y(slf,:); %JFF addition
    Y(sk,:) = Y(sk,:) + F.factors(i).T*Y(rd,:);
    Y(sk,:) = F.factors(i).Ds\Y(sk,:);  %JFF addition
    Y(rd,:) = F.factors(i).L'*Y(rd,:);
    Y(rd,:) = Y(rd,:) + F.factors(i).E'*Y(sk,:);
  end

  % downward sweep
  for i = n:-1:1
    sk = F.factors(i).sk;
    rd = F.factors(i).rd;
    slf = F.factors(i).slf; %JFF addition
    Y(sk,:) = Y(sk,:) + F.factors(i).E*Y(rd,:);
    Y(rd,:) = F.factors(i).L*Y(rd,:);
    Y(sk,:) = F.factors(i).Ds\Y(sk,:);  %JFF addition
    Y(rd,:) = Y(rd,:) + F.factors(i).T'*Y(sk,:);
    Y(slf,:) = F.factors(i).D*Y(slf,:); %JFF addition
  end
end