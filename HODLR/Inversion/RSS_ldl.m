function [F] = RSS_ldl(H,tol,fID)
%HIF Summary of this function goes here
%   Detailed explanation goes here
fprintf(fID,'--------------HIF contruction------------\n');
fprintf(fID,'lvl | DOFS | rd | sk | rd/DOFS | sk/DOFS | time (s)\n')
N = H.N;
n_b = 2^(H.nlvl-1); 
DOFS = true(N,1);
                  
e = cell(n_b,1);
F = struct('slf',e,'rd',e,'sk',e,'L',e,'U',e,'T',e,'C1',e,'C2',e);
F = struct('factors',F,'nlvl',H.nlvl,'N',H.N,'nblocks',0);
n_f = 0;

for lvl=H.nlvl-1:-1:1
    tic
    blk_size = N/2^lvl;
    nDOFS = sum(DOFS);
    for blk=1:2:2^lvl
        %save('check_HIF','lvl','blk','n_f','H')
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
            [L_i,D_i] = ldl(A_ii(rd_i,rd_i)- 1i*diag(imag(diag(A_ii(rd_i,rd_i)))));
            U_i = L_i';
            L_i = L_i*D_i;
            C_i_1 = A_ii(sk_i,rd_i)/U_i;
            C_i_2 = L_i\A_ii(rd_i,sk_i);
            A_ii_sk = A_ii(sk_i,sk_i) - C_i_1*C_i_2;
            DOFS(xi_rem(rd_i))=0;
            
            % Save skeletonization factors in F
            while n_b < n_f+1
              e = cell(n_b,1);
              s  = struct('slf',e,'rd',e,'sk',e,'L',e,'U',e,'T',e,'C1',e,'C2',e);
              F.factors = [F.factors; s];
              n_b = 2*n_b;
            end
            F.factors(n_f).T = T_i;
            F.factors(n_f).sk = xi_rem(sk_i);
            F.factors(n_f).rd = xi_rem(rd_i);
            F.factors(n_f).L = L_i;
            F.factors(n_f).C1 = C_i_1; 
            F.factors(n_f).U = U_i;
            F.factors(n_f).C2 = C_i_2; 
        end
        
        % Do skeletonization in right child if possible
        if ~isempty(rd_j) 
            n_f = n_f + 1;
            A_jj = H.factors(blk+1).A_blkdiag;
            
            A_jj(rd_j,:) = A_jj(rd_j,:) - T_j'*A_jj(sk_j,:);
            A_jj(:,rd_j) = A_jj(:,rd_j) - A_jj(:,sk_j)*T_j;
            [L_j,D_j] = ldl(A_jj(rd_j,rd_j)- 1i*diag(imag(diag(A_jj(rd_j,rd_j)))));
            U_j = L_j';
            L_j = L_j*D_j;
            C_j_1 = A_jj(sk_j,rd_j)/U_j;   
            C_j_2 = L_j\A_jj(rd_j,sk_j); 
            A_jj_sk = A_jj(sk_j,sk_j) - C_j_1*C_j_2;
            DOFS(xj_rem(rd_j))=0;
            F.factors(n_f).T = T_j;
            F.factors(n_f).sk = xj_rem(sk_j);
            F.factors(n_f).rd = xj_rem(rd_j);   
            F.factors(n_f).L = L_j;      
            F.factors(n_f).C1 = C_j_1;
            F.factors(n_f).U = U_j;      
            F.factors(n_f).C2 = C_j_2;
        end
        
        
        
        % merge both child to create parent node
        pos_i = H.factors(blk).pos(2);
        pos_j = H.factors(blk+1).pos(2);
        sk_i = DOFS(xi);
        sk_j = DOFS(xj);
        V_ij = H.factors(blk).V(DOFS(xi),1:pos_i);
        U_ij = H.factors(blk+1).V(DOFS(xj),1:pos_j);
        %save('c_HIF2','A_ii_sk','V_ij','U_ij','A_jj_sk')
        A_ij = [A_ii_sk, V_ij*U_ij'; U_ij*V_ij' A_jj_sk];
        %A_off2 = U_ij*V_ij';
        %save('check_HIF','A_ii','A_jj','A_jj_sk','A_ii_sk','A_ij','U_ij','V_ij','lvl','blk','n_f','sk_j','sk_i','H','A_ir','DOFS','xi','xj','F')
        %if blk==1 && lvl==H.nlvl-1
            %save('check_HIFsym','T_i','T_j','sk_j','sk_i','rd_i','rd_j','A_ij','A_ii','A_jj','L_i','L_j','V_ij','U_ij')
        %end
        blk_merge = (blk-1)/2+1; 
        H.factors(blk_merge).A_blkdiag = A_ij;
        H.factors(blk_merge).slf = [H.factors(blk).slf,H.factors(blk+1).slf];
        if lvl~=1
            %save('c_HIF3','pos_i','pos_j','blk','blk_merge','H')
            H.factors(blk_merge).V = [H.factors(blk).V(:,pos_i+1:end);H.factors(blk+1).V(:,pos_j+1:end)];
            H.factors(blk_merge).U = [H.factors(blk).U(:,pos_i+1:end);H.factors(blk+1).U(:,pos_j+1:end)];
            H.factors(blk_merge).pos = H.factors(blk).pos(2:end)-H.factors(blk).pos(2);
        else
            n_f = n_f+1;
            [L,D] = ldl(A_ij- 1i*diag(imag(diag(A_ij))));
            U = L';
            L = L*D;
            F.factors(n_f).L = L;
            F.factors(n_f).U = U;
            F.factors(n_f).sk = [xi(sk_i),xj(sk_j)];
        end
        
    end
    t = toc;
    nDOFS_rem = sum(DOFS);
    fprintf(fID,'%3d | %6d | %8d | %8d | %8.2f | %8.2f | %10.2e (s)\n', ...
                    lvl,nDOFS,nDOFS-nDOFS_rem,nDOFS_rem,(nDOFS-nDOFS_rem)/nDOFS,nDOFS_rem/nDOFS,toc);
end
F.nblocks = n_f;
fprintf(fID,'-------------End HIF contruction-----------\n');
end

%H_aux = H;

%lvl = nlvl-1;
%n_f = 1;
%for nb=1:2^lvl
%   H_aux(nb).slf = F.factors(nb).slf;
%   H_aux(nb).A_blkdiag = F.factors(nb).A_blkdiag;
%   H_aux(nb).blk = [nb]; 
%end
