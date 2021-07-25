function [ F, HODLR] = HODLR_construction_skew( N, func, tol, fid, occ, mr_1, mr_2)

    %Construct HOLDR representation of the matrix using peeling algorithm. 
    % func handle provides matrix vector multiplication 
    % Low rank matrices are represented using notation A\approx U*V
    % structure factors contains the information of the matrix. For each block
    % at the leave level, slf contains the indices for DOFs, U and V the
    % corresponding lowrank matrices at each level 2:L for indices slf, and pos
    % indicates which columns of U and V correspond to which level.
    
    nlvl = ceil(max(0,log2(N/occ)))+1;
    
    mn = 2^(nlvl-1);                  
    e = cell(mn,1);
    F = struct('slf',e,'U',e,'V',e,'pos',0,'A_blkdiag',e);
    F = struct('factors',F,'nlvl',nlvl,'N',N);
    %Ftemp = struct('S',e,'lvl',e,'blk',e);
    
    
    HODLR = struct('x',e,'y',e,'U',e,'V',e,'S',e,'nb',cell(1,1),'nd',cell(1,1),'lvl',e);
             
    nf = mn;
    nb = 0;
    c = 2;
    leaf_size = N/2^(nlvl-1);
    
    nS = mn;
    
    fprintf(fid,'------------------HODLR construction-----------------\n');
    fprintf(fid,'level |  rank  | maxrank | time (s)\n');
    for lvl=1:nlvl-1
       tic
       if lvl==1
          max_rank = min(mr_1,N/2^lvl-c);
       else
          max_rank = min(mr_2,N/2^lvl-c);
       end
       R_1 = zeros(N,max_rank+c);
       R_2 = zeros(N,max_rank+c);
       % Generate random matrices
       for blk=1:2:2^lvl
           blk_size = N/2^lvl;
           R_1((blk-1)*blk_size+1:blk*blk_size,:)=random('normal',0,1,[blk_size,max_rank+c]);
           R_2(blk*blk_size+1:(blk+1)*blk_size,:)=random('normal',0,1,[blk_size,max_rank+c]);
       end
    
       %Apply A to random matrices and substract H-matrix constructed so far
       if lvl>1
            A_1 = func(R_1)-HODLR_apply_skew(HODLR,R_1);
            A_2 = func(R_2)-HODLR_apply_skew(HODLR,R_2);
       else
           A_1 = func(R_1);
           A_2 = func(R_2);
       end
       Q_row = zeros(N,max_rank+c);
       for blk=1:2:2^lvl
           [Q_row((blk-1)*blk_size+1:blk*blk_size,:),~,~] = qr(A_2((blk-1)*blk_size+1:blk*blk_size,:),0);
       end
       if lvl>1
            AQ_row = func(Q_row)-HODLR_apply_skew(HODLR,Q_row);
       else
           AQ_row = func(Q_row);
       end
       
       nS = nS - 2^(lvl-1); 
       for blk=1:2:2^lvl
           nb = nb + 1;
           AR1 = A_1(blk*blk_size+1:(blk+1)*blk_size,:);
           [Q_col,~,~] = qr(AR1,0);
           AQr = AQ_row(blk*blk_size+1:(blk+1)*blk_size,:);
           QcAQr = Q_col'*AQr;
           
           U=[];V=[];
           [U,S,V] = RSVD(QcAQr,tol,max_rank);
           if ~isempty(U)
              U = Q_col*U;
              V = Q_row((blk-1)*blk_size+1:blk*blk_size,:)*V;
              
           end
          
           nf = nf + 1;
           while mn < nf
              e = cell(mn,1);
              s = struct('x',e,'y',e,'U',e,'V',e,'S',e,'nb',e,'nd',cell(1,1),'lvl',e);
              HODLR = [HODLR; s];
              mn = 2*mn;
           end
           
           HODLR(nb).x = blk*blk_size+1:(blk+1)*blk_size;
           HODLR(nb).y = (blk-1)*blk_size+1:blk*blk_size;
           HODLR(nb).U = U;
           HODLR(nb).V = V;
           HODLR(nb).S = S;
           HODLR(nb).lvl = lvl;
           HODLR(1).nb = nb;
           
           nS = nS+1;
           %Ftemp(nS).S = S;
           %Ftemp(nS).blk = blk;
           %Ftemp(nS).lvl = lvl;
    
           for sblk=1:2^(nlvl-lvl-1)
              idx_blk = (blk-1)*2^(nlvl-lvl-1)+sblk;
              F.factors(idx_blk).V = [F.factors(idx_blk).V,V((sblk-1)*leaf_size+1:sblk*leaf_size,:)];
              F.factors(idx_blk).U = [F.factors(idx_blk).U,U((sblk-1)*leaf_size+1:sblk*leaf_size,:)*S];          
              F.factors(idx_blk).pos = [F.factors(idx_blk).pos,F.factors(idx_blk).pos(end)+size(U,2)];
    
              F.factors(idx_blk+2^(nlvl-lvl-1)).V = [F.factors(idx_blk+2^(nlvl-lvl-1)).V,U((sblk-1)*leaf_size+1:sblk*leaf_size,:)];
              F.factors(idx_blk+2^(nlvl-lvl-1)).U = [F.factors(idx_blk+2^(nlvl-lvl-1)).U,V((sblk-1)*leaf_size+1:sblk*leaf_size,:)*S];
              F.factors(idx_blk+2^(nlvl-lvl-1)).pos = [F.factors(idx_blk+2^(nlvl-lvl-1)).pos,F.factors(idx_blk+2^(nlvl-lvl-1)).pos(end)+size(U,2)];
           end
       end
       nS = nS - 2^(lvl-1); 
       fprintf(fid,'%3d | %6d | %4d | %10.2e (s)\n',lvl,size(U,2),max_rank,toc);   
    end
    %F.Sm(1:length(Ftemp)-nS) = Ftemp(nS+1:length(Ftemp));
    
    nd = nb;
    blk_size = N/2^lvl;
    R_1 = repmat(eye(blk_size),[2^lvl,1]);
    
    A_1 = func(R_1)-HODLR_apply_skew(HODLR,R_1);
    for blk=1:2^(lvl)
        nd = nd+1;
        HODLR(nd).U = (A_1((blk-1)*blk_size+1:blk*blk_size,:)-A_1((blk-1)*blk_size+1:blk*blk_size,:)')/2;
        % HODLR(nd).U = A_1((blk-1)*blk_size+1:blk*blk_size,:);
        HODLR(nd).x = (blk-1)*blk_size+1:blk*blk_size;
        HODLR(nd).lvl = lvl;
        F.factors(blk).A_blkdiag = (A_1((blk-1)*blk_size+1:blk*blk_size,:)-A_1((blk-1)*blk_size+1:blk*blk_size,:)')/2;
        % F.factors(blk).A_blkdiag = A_1((blk-1)*blk_size+1:blk*blk_size,:);
        F.factors(blk).slf = (blk-1)*blk_size+1:blk*blk_size;
        pos = flip(F.factors(blk).pos);
        pos_=zeros(size(pos));
        for i=2:length(pos)
            pos_(i) = pos(i-1)-pos(i)+pos_(i-1);
        end
        F.factors(blk).pos = pos_;
        F.factors(blk).V = flip(F.factors(blk).V,2);
        F.factors(blk).U = flip(F.factors(blk).U,2);
    end
    HODLR(1).nd = nd;
    end