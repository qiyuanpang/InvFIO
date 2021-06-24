% HIFIE2   Constructs hierarchical interpolative factorization for matrix
%          represeneted in Hmatrix format.
%
%    This code is extracted from https://github.com/klho/FLAM and modified 
%    for Hmatrix format.
%
%
%    [F,t] = HIF_Hmatrix(H_matrix,X,rank_or_tol,t,opts,fileID)
%    produces a factorization F of the interaction matrix A on the points X
%    represented in the H-matrix format stored in H_matrix. It uses local 
%    precision parameter RANK_OR_TOL, and proxy function PXYFUN to capture
%    the far field. This is a function of the form
%
%    Parameter opts allows to define:
%
%      - SKIP: skip the dimension reductions on the first SKIP levels (default:
%              SKIP = 0).
%
%      - SYMM: assume that the matrix is unsymmetric if SYMM = 'N', (complex-)
%              symmetric if SYMM = 'S', Hermitian if SYMM = 'H', and Hermitian
%              positive definite if SYMM = 'P' (default: SYMM = 'N'). If
%              SYMM = 'N' or 'S', then local factors are computed using the LU
%              decomposition; if SYMM = 'H', the LDL decomposition; and if
%              SYMM = 'P', the Cholesky decomposition.
%
%      - VERB: display status of the code if VERB = 1 (default: VERB = 0).
%
%    References:
%
%      E. Corona, P.-G. Martinsson, D. Zorin. An O(N) direct solver for
%        integral equations on the plane. Appl. Comput. Harmon. Anal. 38 (2):
%        284-317, 2015.
%
%      K.L. Ho, L. Ying. Hierarchical interpolative factorization for elliptic
%        operators: integral equations. Comm. Pure Appl. Math. 69 (7):
%        1314-1353, 2016.


function [F,t] = HIF_Hmatrix(H_matrix,x,rank_or_tol,t,opts,fileID)
  start = tic;

  % set default parameters
  if nargin < 5
    pxyfun = [];
  end
  if nargin < 6
    opts = [];
  end
  if ~isfield(opts,'ext')
    opts.ext = [];
  end
  if ~isfield(opts,'lvlmax')
    opts.lvlmax = Inf;
  end
  if ~isfield(opts,'skip')
    opts.skip = 0;
  end
  if ~isfield(opts,'symm')
    opts.symm = 'n';
  end
  if ~isfield(opts,'verb')
    opts.verb = 0;
  end

  % check inputs
  assert(opts.skip >= 0,'FLAM:hifie2:negativeSkip', ...
         'Skip parameter must be nonnegative.')
  assert(strcmpi(opts.symm,'n') || strcmpi(opts.symm,'s') || ...
         strcmpi(opts.symm,'h') || strcmpi(opts.symm,'p'), ...
         'FLAM:hifie2:invalidSymm', ...
         'Symmetry parameter must be one of ''N'', ''S'', ''H'', or ''P''.')

  N = size(x,2);

% pull up skeletons for leaves
  for i = t.lvp(t.nlvl)+1:t.lvp(t.nlvl+1)
    node_size = size(t.nodes(i).xi,2);
    t.nodes(i).sk_DOFS = 1:node_size;
  end
  for i = t.lvp(1)+1:t.lvp(t.nlvl+1)
      t.nodes(i).update = false;
  end

  % count nonempty boxes at each level
  pblk = zeros(t.nlvl+1,1);
  for lvl = 1:t.nlvl
    pblk(lvl+1) = pblk(lvl);
    for i = t.lvp(lvl)+1:t.lvp(lvl+1)
      if ~isempty(t.nodes(i).xi)
        pblk(lvl+1) = pblk(lvl+1) + 1;
      end
    end
  end

  % initialize
  mn = t.lvp(end);
  e = cell(mn,1);
  F = struct('sk',e,'rd',e,'T',e,'E',e,'F',e,'L',e,'U',e);
  F = struct('N',N,'nlvl',0,'lvp',zeros(1,2*t.nlvl+1),'factors',F,'symm', ...
             opts.symm);
  nlvl = 0;
  n = 0;
  rem_ = true(N,1);
  mnz = 128;
  M = sparse(N,N);
  
  nz_B = 0;
  mnz_B = 128;
  B = sparse(N,N);
  I_B = zeros(mnz,1);
  J_B = zeros(mnz,1);
  S_B = zeros(mnz,1);
  for node = t.lvp(t.nlvl)+1:t.lvp(t.nlvl+1)
      for node_nbr = [node, t.nodes(node).nbor]
          S_B_ = H_matrix{t.nlvl+1}(node-t.lvp(t.nlvl), node_nbr-t.lvp(t.nlvl)).M; 
          node_DOFs = t.nodes(node).xi;
          node_nbr_DOFs = t.nodes(node_nbr).xi;
          [I_B_,J_B_] = ndgrid(node_DOFs, node_nbr_DOFs);
          m = length(I_B_(:));
          while mnz_B < nz_B + m
            e = zeros(mnz_B,1);
            I_B = [I_B; e];
            J_B = [J_B; e];
            S_B = [S_B; e];
            mnz_B = 2*mnz_B;
          end
          I_B(nz_B+1:nz_B+m) = I_B_(:);
          J_B(nz_B+1:nz_B+m) = J_B_(:);
          S_B(nz_B+1:nz_B+m) = S_B_(:);
          nz_B = nz_B + m;
      end
  end  
  B = sparse(I_B(1:nz_B),J_B(1:nz_B),S_B(1:nz_B),N,N);

  I = zeros(mnz,1);
  J = zeros(mnz,1);
  S = zeros(mnz,1);
  P = zeros(N,1);
  fprintf(fileID,'nlvl: %6d \n', t.nlvl);
  % loop over tree levels
  
  lvl_thrshold = max(t.nlvl-1,3);
  for lvl = t.nlvl:-1:1
    l = t.lrt/2^(lvl - 1);
    nbox = t.lvp(lvl+1) - t.lvp(lvl);
    
    % pull up skeletons from children 
    for i = t.lvp(lvl)+1:t.lvp(lvl+1)
       chld_size = size(t.nodes(i).xi,2)/4;
       chldrn = t.nodes(i).chld;
      for idx_child = 1:length(chldrn)
        child = chldrn(idx_child);
        t.nodes(i).sk_DOFS = [t.nodes(i).sk_DOFS, t.nodes(child).sk_DOFS...
            + chld_size*(idx_child-1)];
      end
    end

    % loop over dimensions
    for d = [2,1]
      tic

      % dimension reduction
      if d < 2

        % continue if in skip stage
        if lvl > t.nlvl - opts.skip
          continue
        end
        
        % generate edge centers
        ctr = zeros(4*nbox,2);
        box2ctr = cell(nbox,1);
        for i = t.lvp(lvl)+1:t.lvp(lvl+1)
          j = i - t.lvp(lvl);
          idx = 4*(j-1)+1:4*j;
          off = [0 -1; -1  0; 0 1; 1 0];
          ctr(idx,:) = bsxfun(@plus,t.nodes(i).ctr,0.5*l*off);
          box2ctr{j} = idx;
        end

        % find unique shared centers
        idx = bsxfun(@minus,ctr,t.nodes(1).ctr);
        idx = round(2*idx/l);
        max_idx = max(max(idx));
        idx(idx==max_idx) = -max_idx;

        % initialize
        nb = size(ctr,1);
        e = cell(nb,1);
        blocks_ = struct('ctr',e,'xi',e,'prnt',e,'sk',e);
        e = cell(nb/2,1);
        blocks = struct('ctr',e,'xi',e,'prnt_1',e,'prnt_2',e,'sk_1',e,'sk_2',e);
        for i = 1:nb
          blocks_(i).ctr = ctr(i,:);
        end
        
        % sort points by centers
        for box = 1:nbox
          sk_DOFS = t.nodes(t.lvp(lvl)+box).sk_DOFS;
          xi = [t.nodes(t.lvp(lvl)+box).xi(t.nodes(t.lvp(lvl)+box).sk_DOFS)];
          i = box2ctr{box};
          dx = bsxfun(@minus,x(1,xi),ctr(i,1));
          dy = bsxfun(@minus,x(2,xi),ctr(i,2));
          dist = sqrt(dx.^2 + dy.^2);
          near = bsxfun(@eq,dist,min(dist,[],1));
          for i = 1:length(xi)
            j = find(near(:,i),1);
            blocks_(box2ctr{box}(j)).sk = ...
                [blocks_(box2ctr{box}(j)).sk, sk_DOFS(i)];
            blocks_(box2ctr{box}(j)).xi = ...
                [blocks_(box2ctr{box}(j)).xi, xi(i)];
            blocks_(box2ctr{box}(j)).prnt = t.lvp(lvl)+box;
          end
        end
        
        [~,i,j] = unique(idx,'rows');
        for edge = 1:size(idx,1)
            unique_edge = j(edge);
            blocks(unique_edge).ctr = blocks_(edge).ctr;
            blocks(unique_edge).xi = [blocks(unique_edge).xi, blocks_(edge).xi];
            if isempty(blocks(unique_edge).sk_1)
                blocks(unique_edge).sk_1 = blocks_(edge).sk;
                blocks(unique_edge).prnt_1 = blocks_(edge).prnt;
            else
                blocks(unique_edge).sk_2 = blocks_(edge).sk;
                blocks(unique_edge).prnt_2 = blocks_(edge).prnt;
            end               
        end
      end

      % initialize
      nlvl = nlvl + 1;
      if d == 2
        nb = t.lvp(lvl+1) - t.lvp(lvl);
      else
        nb = length(blocks);
      end
      nrem1 = sum(rem_);
      nblk = pblk(lvl) + nb;
      nz = 0;
      nz_B = 0;

      % loop over blocks
      for i = 1:nb
        if d == 2
          j = t.lvp(lvl) + i;
          blk = t.nodes(j);
          nbr = [t.nodes(blk.nbor).xi];
          sk_DOFS = blk.sk_DOFS;
          slf = blk.xi(sk_DOFS);
          nslf = length(slf);
          sslf = sort(slf);
        else
          blk = blocks(i);
          nbor = [t.nodes(blk.prnt_1).nbor,t.nodes(blk.prnt_2).nbor];
          nbr = unique([t.nodes(nbor).xi]);
          nbr(ismember(nbr,blk.xi))=[];
          nbr = nbr(rem_(nbr)>0);
          sk_DOFS_1 = blk.sk_1;
          sk_DOFS_2 = blk.sk_2;
          slf = blk.xi;
          nslf = length(slf);
          sslf = sort(slf);
        end

        % add neighbors with modified interactions
        [mod,~] = find(M(:,slf));
        mod = unique(mod);
        mod = mod(~ismembc(mod,sslf));
        nbr = unique([nbr(:); mod(:)]);
        nnbr = length(nbr);
        snbr = sort(nbr);

        % compute interaction matrix
        if d==2
            % Obtain near field
            K = full(B(nbr,slf)) + spget('nbr','slf');
            
            % Add far-field low-rank bases
            if lvl>lvl_thrshold
                K2 = Hget(sk_DOFS,lvl,t,i);
                K = [K; K2];
            end
            
        else
            % Obtain near field
            K = full(B(nbr,slf)) + spget('nbr','slf');%full(M(nbr,slf)); %spget('nbr','slf');
            
           % Add far-field low-rank bases
            if lvl>lvl_thrshold
                K2 = Hget(sk_DOFS_1,lvl,t,blk.prnt_1 - t.lvp(lvl));
                K3 = Hget(sk_DOFS_2,lvl,t,blk.prnt_2 - t.lvp(lvl));
                K = [K; blkdiag(K2,K3)];
            end
        end
        
        
        % skeletonize
        [sk,rd,T] = id(K,rank_or_tol,1);
        % restrict to skeletons
        if d == 2
          t.nodes(j).sk_DOFS = sort(sk_DOFS(sk)); 
        else  
          sk_1 = sk(sk<=length(sk_DOFS_1));
          sk_2 = sk(sk>length(sk_DOFS_1)) - length(sk_DOFS_1);
          if t.nodes(blk.prnt_1).update
              t.nodes(blk.prnt_1).sk_DOFS = sort([t.nodes(blk.prnt_1).sk_DOFS, sk_DOFS_1(sk_1)]);
          else
              t.nodes(blk.prnt_1).sk_DOFS = sort(sk_DOFS_1(sk_1));
              t.nodes(blk.prnt_1).update = true;
          end
          if t.nodes(blk.prnt_2).update
              t.nodes(blk.prnt_2).sk_DOFS = sort([t.nodes(blk.prnt_2).sk_DOFS, sk_DOFS_2(sk_2)]);
          else
              t.nodes(blk.prnt_2).sk_DOFS = sort(sk_DOFS_2(sk_2));
              t.nodes(blk.prnt_2).update = true;
          end
        end

        % move on if no compression
        if isempty(rd)
          continue
        end
        rem_(slf(rd)) = 0;

        % compute factors
        K = full(B(slf,slf)) + spget('slf','slf');

        if strcmpi(opts.symm,'s')
          K(rd,:) = K(rd,:) - T.'*K(sk,:);
        else
          K(rd,:) = K(rd,:) - T'*K(sk,:);
        end
        K(:,rd) = K(:,rd) - K(:,sk)*T;
        if strcmpi(opts.symm,'n') || strcmpi(opts.symm,'s')
          [L,U] = lu(K(rd,rd));
          E = K(sk,rd)/U;
          G = L\K(rd,sk);
        elseif strcmpi(opts.symm,'h')
          B_rd = K(rd,rd);
          B_rd = B_rd - diag(diag(B_rd)) + diag(real(diag(B_rd)));
          [L,U] = ldl(B_rd);
          E = (K(sk,rd)/L')/U;
          G = [];
        elseif strcmpi(opts.symm,'p')
          L = chol(K(rd,rd),'lower');
          U = [];
          E = K(sk,rd)/L';
          G = [];
        end

        % update self-interaction
        if strcmpi(opts.symm,'n') || strcmpi(opts.symm,'s')
          S_ = -E*G;
        elseif strcmpi(opts.symm,'h')
          S_ = -E*U*E';
        elseif strcmpi(opts.symm,'p')
          S_ = -E*E';
        end
        [I_,J_] = ndgrid(slf(sk));
        m = length(sk)^2;
        while mnz < nz + m
          e = zeros(mnz,1);
          I = [I; e];
          J = [J; e];
          S = [S; e];
          mnz = 2*mnz;
        end
        I(nz+1:nz+m) = I_(:);
        J(nz+1:nz+m) = J_(:);
        S(nz+1:nz+m) = S_(:);
        nz = nz + m;

        % store matrix factors
        n = n + 1;
        while mn < n
          e = cell(mn,1);
          s = struct('sk',e,'rd',e,'T',e,'E',e,'F',e,'L',e,'U',e);
          F.factors = [F.factors; s];
          mn = 2*mn;
        end
        F.factors(n).sk = slf(sk);
        F.factors(n).rd = slf(rd);
        F.factors(n).T = T;
        F.factors(n).E = E;
        F.factors(n).F = G;
        F.factors(n).L = L;
        F.factors(n).U = U;
      end
      F.lvp(nlvl+1) = n;

      % update modified entries
      [I_,J_,S_] = find(M);
      idx = rem_(I_) & rem_(J_);
      I_ = I_(idx);
      J_ = J_(idx);
      S_ = S_(idx);
      m = length(S_);
      while mnz < nz + m
        e = zeros(mnz,1);
        I = [I; e];
        J = [J; e];
        S = [S; e];
        mnz = 2*mnz;
      end
      I(nz+1:nz+m) = I_;
      J(nz+1:nz+m) = J_;
      S(nz+1:nz+m) = S_;
      nz = nz + m;
      M = sparse(I(1:nz),J(1:nz),S(1:nz),N,N);
      
      % update B matrix
      
      if d==2
           % loop over blocks
          for i = 1:nb
            %Store near-field interaction of skeletons
            if lvl>lvl_thrshold
                node = i + t.lvp(lvl);
                for node_nbr = t.nodes(node).inter
                  j = node_nbr - t.lvp(lvl);
                  if node_nbr <= node
                      continue
                  end
                  xi_sk = t.nodes(node).sk_DOFS;
                  xj_sk = t.nodes(node_nbr).sk_DOFS;
                  node_DOFs = t.nodes(node).xi(xi_sk);
                  node_nbr_DOFs = t.nodes(node_nbr).xi(xj_sk);
                  if ~isempty(H_matrix{lvl,1}(i,j).U)
                    S_B_ = H_matrix{lvl,1}(i,j).U(xi_sk,:)...
                     *H_matrix{lvl,1}(i,j).M*H_matrix{lvl,1}(j,i).U(xj_sk,:)';
                  else
                      continue
                  end
                  [I_B_,J_B_] = ndgrid(node_DOFs, node_nbr_DOFs);
                  m = length(I_B_(:));
                  while mnz_B < nz_B + 2*m
                    e = zeros(mnz_B,1);
                    I_B = [I_B; e];
                    J_B = [J_B; e];
                    S_B = [S_B; e];
                    mnz_B = 2*mnz_B;
                  end
                  I_B(nz_B+1:nz_B+m) = I_B_(:);
                  J_B(nz_B+1:nz_B+m) = J_B_(:);
                  S_B(nz_B+1:nz_B+m) = S_B_(:);
                  S_B_ = S_B_'; I_B_ = I_B_'; J_B_ = J_B_';
                  I_B(nz_B+m+1:nz_B+2*m) = J_B_(:);
                  J_B(nz_B+m+1:nz_B+2*m) = I_B_(:);
                  S_B(nz_B+m+1:nz_B+2*m) = S_B_(:);
                  nz_B = nz_B + 2*m;
                end
            end 
          end

          if lvl>lvl_thrshold
              [I_B_,J_B_,S_B_] = find(B);
              idx = rem_(I_B_) & rem_(J_B_);
              I_B_ = I_B_(idx);
              J_B_ = J_B_(idx);
              S_B_ = S_B_(idx);
              m = length(S_B_);
              while mnz_B < nz_B + m
                e = zeros(mnz_B,1);
                I_B = [I_B; e];
                J_B = [J_B; e];
                S_B = [S_B; e];
                mnz_B = 2*mnz_B;
              end
              I_B(nz_B+1:nz_B+m) = I_B_;
              J_B(nz_B+1:nz_B+m) = J_B_;
              S_B(nz_B+1:nz_B+m) = S_B_;
              nz_B = nz_B + m;

              B = sparse(I_B(1:nz_B),J_B(1:nz_B),S_B(1:nz_B),N,N);
          end
      end

      % print summary
      if opts.verb
        nrem2 = sum(rem_);
        fprintf(fileID,'%3d-%1d | %6d | %8d | %8d | %8.2f | %8.2f | %10.2e (s)\n', ...
                lvl,d,nb,nrem1,nrem2,nrem1/nb,nrem2/nb,toc);
      end
      if nblk == 1
        break
      end
    end
  end

  % finish
  F.nlvl = nlvl;
  F.lvp = F.lvp(1:nlvl+1);
  F.factors = F.factors(1:n);
  if opts.verb
    fprintf(fileID,['-'*ones(1,80) '\n']);
    wtime = toc(start);
    fprintf (fileID, 'Elapsed time is %f seconds.\n', wtime);
  end

  % sparse matrix access function (native MATLAB is slow for large matrices)
  function A = spget(Ityp,Jtyp)
    if strcmpi(Ityp,'slf')
      I_ = slf;
      m_ = nslf;
      I_sort = sslf;
    elseif strcmpi(Ityp,'nbr')
      I_ = nbr;
      m_ = nnbr;
      I_sort = snbr;
    end
    if strcmpi(Jtyp,'slf')
      J_ = slf;
      n_ = nslf;
    elseif strcmpi(Jtyp,'nbr')
      J_ = nbr;
      n_ = nnbr;
    end
    P(I_) = 1:m_;
    A = zeros(m_,n_);
    [I_,J_,S_] = find(M(:,J_));
    idx = ismembc(I_,I_sort);
    I_ = I_(idx);
    J_ = J_(idx);
    S_ = S_(idx);
    A(P(I_) + (J_ - 1)*m_) = S_;
  end

%  matrix access function for H-matrix 
  function B = Hget(sk_nodes,lvl,t,node)
    B = [];

    for level = lvl:-1:lvl_thrshold+1
        node_level = ceil(node/4^(lvl-level));
        local_node = rem(node,4^(lvl-level));
        if local_node == 0
            local_node = 4^(lvl-level);
        end
        blk_size = size(t.nodes(node_level + t.lvp(level)).xi,2)/(4^(lvl-level));
        idx_nodes = blk_size*(local_node-1) + sk_nodes;
        for nbr_node=t.nodes(node_level + t.lvp(level)).inter
           if node_level+t.lvp(level)<nbr_node
            B = [B; H_matrix{level,1}(node_level,nbr_node - t.lvp(level)).M'*H_matrix{level,1}(node_level,nbr_node - t.lvp(level)).U(idx_nodes,:)'];
           else
            B = [B; H_matrix{level,1}(nbr_node - t.lvp(level),node_level).M*H_matrix{level,1}(node_level,nbr_node - t.lvp(level)).U(idx_nodes,:)'];   
           end
        end   
    end
  end
end
