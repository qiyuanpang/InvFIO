function [sk,rd,T] = ID(A,rank_or_tol,srand)
% ID   Interpolative decomposition.
%
%    [SK,RD,T] = ID(A,tol) produces a rank-K approximation of A via the skeleton
%    and redundant indices SK and RD, respectively, and an interpolation matrix
%    T such that A(:,RD) = A(:,SK)*T + E, where LENGTH(SK) = K and NORM(E)
%    is of order TOL*NORM(A).
%
%    Obtain from FLAM: https://github.com/klho/FLAM

  % set default parameters
  if nargin < 3 || isempty(srand)
    srand = 1;
  end

  % check inputs
  assert(rank_or_tol >= 0,'FLAM:id:negativeRankOrTol', ...
         'Rank or tolerance must be nonnegative.')

  % initialize
  [m,n] = size(A);

  % return if matrix is empty
  if isempty(A)
    sk = [];
    rd = 1:n;
    T = zeros(0,n);
    return
  end

  % sample against Gaussian matrix if too rectangular
  if srand && m > 2*n
    A = randn(n+16,m)*A;
  end

  % compute ID
  [~,R,E] = qr(A,0);
  if rank_or_tol < 1
    k = sum(abs(diag(R)) > abs(R(1))*rank_or_tol);
  else
    k = min(rank_or_tol,n);
  end
  [sk,psk] = sort(E(1:k));
  [rd,psd] = sort(E(k+1:end));
  T = R(1:k,1:k)\R(1:k,k+1:end);
  T = T(psk,psd);
end
