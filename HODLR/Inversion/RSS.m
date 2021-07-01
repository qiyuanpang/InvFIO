function [F] = RSS(H,tol,fID)

fprintf(fID,'--------------RSS contruction 1D------------\n');
fprintf(fID,'lvl | DOFS | rd | sk | rd/DOFS | sk/DOFS | time (s)\n');
N = H.N;
n_b = 2^(H.nlvl-1); 
DOFS = true(N,1);
                  
e = cell(n_b,1);
F = struct('slf',e,'rd',e,'sk',e,'L',e,'T',e,'C',e);
F = struct('factors',F,'nlvl',H.nlvl,'N',H.N,'nblocks',0);
n_f = 0;

for lvl=H.nlvl-1:-1:1
    tic
    blk_size = N/2^lvl;
    nDOFS = sum(DOFS);
    for blk=1:2:2^lvl
        xi = (blk-1)*blk_size+1:blk*blk_size;
        xj = blk*blk_size+1:(blk+1)*blk_size;

        % Find skeletons
        xi_rem = xi(DOFS(xi));
        xj_rem = xj(DOFS(xj));
        A_ir = H.factors(blk).V(DOFS(xi),:)';  
        A_jr = H.factors(blk+1).V(DOFS(xj),:)';
        [sk_i,rd_i,T_i] = ID(A_ir,tol);
        [sk_j,rd_j,T_j] = ID(A_jr,tol);
        A_ii_sk = H.factors(blk).A_blkdiag;
        A_jj_sk = H.factors(blk+1).A_blkdiag;
        
        % Do skeletonization in left child if possible
        if ~isempty(rd_i)
            n_f = n_f + 1;
            A_ii = H.factors(blk).A_blkdiag;
            
            

            A_ii(rd_i,:) = A_ii(rd_i,:) - T_i'*A_ii(sk_i,:);
            A_ii(:,rd_i) = A_ii(:,rd_i) - A_ii(:,sk_i)*T_i;
            % disp(real(eig(A_ii(rd_i,rd_i)- 1i*diag(imag(diag(A_ii(rd_i,rd_i))))))')
            L_i = chol(A_ii(rd_i,rd_i)- 1i*diag(imag(diag(A_ii(rd_i,rd_i)))),'lower');
            % L_i = chol(A_ii(rd_i,rd_i),'lower');
            C_i = A_ii(sk_i,rd_i)/L_i';
            A_ii_sk = A_ii(sk_i,sk_i) - C_i*C_i';
            DOFS(xi_rem(rd_i))=0;
            
            % Save skeletonization factors in F
            while n_b < n_f+1
              e = cell(n_b,1);
              s  = struct('slf',e,'rd',e,'sk',e,'L',e,'T',e,'C',e);
              F.factors = [F.factors; s];
              n_b = 2*n_b;
            end
            F.factors(n_f).T = T_i;
            F.factors(n_f).sk = xi_rem(sk_i);
            F.factors(n_f).rd = xi_rem(rd_i);
            F.factors(n_f).L = L_i;
            F.factors(n_f).C = C_i; 
        end
        
        % Do skeletonization in right child if possible
        if ~isempty(rd_j) 
            n_f = n_f + 1;
            A_jj = H.factors(blk+1).A_blkdiag;
            
            A_jj(rd_j,:) = A_jj(rd_j,:) - T_j'*A_jj(sk_j,:);
            A_jj(:,rd_j) = A_jj(:,rd_j) - A_jj(:,sk_j)*T_j;
            % disp(real(eig(A_jj(rd_j,rd_j)-1i*diag(imag(diag(A_jj(rd_j,rd_j))))))')
            L_j = chol(A_jj(rd_j,rd_j)-1i*diag(imag(diag(A_jj(rd_j,rd_j)))),'lower');
            % L_j = chol(A_jj(rd_j,rd_j),'lower'); 

            C_j = A_jj(sk_j,rd_j)/L_j';        
            A_jj_sk = A_jj(sk_j,sk_j) - C_j*C_j';
            DOFS(xj_rem(rd_j))=0;
            F.factors(n_f).T = T_j;
            F.factors(n_f).sk = xj_rem(sk_j);
            F.factors(n_f).rd = xj_rem(rd_j);   
            F.factors(n_f).L = L_j;      
            F.factors(n_f).C = C_j;
        end
        
        
        
        % merge both child to create parent node
        pos_i = H.factors(blk).pos(2);
        pos_j = H.factors(blk+1).pos(2);
        sk_i = DOFS(xi);
        sk_j = DOFS(xj);
        V_ij = H.factors(blk).V(DOFS(xi),1:pos_i);
        U_ij = H.factors(blk).U(DOFS(xj),1:pos_j);
        A_ij = [A_ii_sk, V_ij*U_ij'; U_ij*V_ij' A_jj_sk];
        
        blk_merge = (blk-1)/2+1; 
        H.factors(blk_merge).A_blkdiag = A_ij;
        H.factors(blk_merge).slf = [H.factors(blk).slf,H.factors(blk+1).slf];
        if lvl~=1
            H.factors(blk_merge).V = [H.factors(blk).V(:,pos_i+1:end);H.factors(blk+1).V(:,pos_j+1:end)];
            H.factors(blk_merge).U = [H.factors(blk).U(:,pos_i+1:end);H.factors(blk+1).U(:,pos_j+1:end)];
            H.factors(blk_merge).pos = H.factors(blk).pos(2:end)-H.factors(blk).pos(2);
        else
            n_f = n_f+1;
            F.factors(n_f).L = chol(A_ij- 1i*diag(imag(diag(A_ij))),'lower');
            F.factors(n_f).sk = [xi(sk_i),xj(sk_j)];
        end
        
    end
    t = toc;
    nDOFS_rem = sum(DOFS);
    fprintf(fID,'%3d | %6d | %8d | %8d | %8.2f | %8.2f | %10.2e (s)\n', ...
                    lvl,nDOFS,nDOFS-nDOFS_rem,nDOFS_rem,(nDOFS-nDOFS_rem)/nDOFS,nDOFS_rem/nDOFS,toc);
end
F.nblocks = n_f;
fprintf(fID,'-------------End RSS contruction 1D -----------\n');
end


