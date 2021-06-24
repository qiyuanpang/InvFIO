function [ T ] = quadtree_strong( N, occ )
%QUADTREE_STRONG Summary of this function goes here
%   Detailed explanation goes here
nlvl = ceil(max(0,log2(N/occ)))+1;
mn = 2^(nlvl-1);                

s = struct('lvl',1,'pos',[],'DOFS',[],'prnt',[],'chld',[],'sbls', [], 'sbl_l',[],...
    'sbl_r',[],'sbl_d',[],'sbl_u',[],'DOFS_l',[],'DOFS_r',[],'DOFS_d',[],'DOFS_u',...
    [],'sbl_n_l',[],'sbl_n_r',[],'sbl_n_d',[],'sbl_n_u',[]);
T = struct('nlvl',nlvl,'nodes',s,'n_blocks',1,'lvp',[]);


T.nodes(2).sbl_r = [3]; T.nodes(2).sbl_d = [4]; T.nodes(2).sbl_l = [3]; T.nodes(2).sbl_u = [4];
T.nodes(3).sbl_l = [2]; T.nodes(3).sbl_d = [5]; T.nodes(3).sbl_r = [2]; T.nodes(3).sbl_u = [5];
T.nodes(4).sbl_u = [2]; T.nodes(4).sbl_r = [5]; T.nodes(4).sbl_d = [2]; T.nodes(4).sbl_l = [5];
T.nodes(5).sbl_l = [4]; T.nodes(5).sbl_u = [3]; T.nodes(5).sbl_r = [4]; T.nodes(5).sbl_d = [3];
DOFS_node = reshape(1:N^2,[N,N])';
T.nodes(1).DOFS = 1:N^2;
for j=1:2
    for k=1:2
        i = 2*(j-1)+k;
        Idx_x = (j-1)*N/2+1:N/2*j;
        Idx_y = (k-1)*N/2+1:N/2*k;
        T.nodes(i+1).DOFS = reshape(DOFS_node(Idx_x,Idx_y)',[N^2/4,1]);
    end
end

n_idx = 0;
T.lvp = [0];

n_lvl = 0;
for lvl=0:nlvl-1
    n_lvl = n_lvl+4^lvl;
    T.lvp = [T.lvp n_lvl];
    for nk=1:4^lvl
        n_size = N/2^lvl;
        node = nk+n_idx;
        chldrn = [n_lvl+(nk-1)*4+1:n_lvl+(nk)*4];
        T.nodes(node).chld = chldrn;
        T.nodes(node).lvl = lvl;
        for j=1:2
            for k=1:2
            i = 2*(j-1)+k;
            DOFS_node = reshape(T.nodes(node).DOFS,[n_size,n_size])';
            Idx_x = (j-1)*n_size/2+1:n_size/2*j;
            Idx_y = (k-1)*n_size/2+1:n_size/2*k;
            DOFS_child = DOFS_node(Idx_x,Idx_y)';
            
            T.nodes(chldrn(i)).DOFS = reshape(DOFS_child,[n_size^2/4,1]);
            T.nodes(chldrn(i)).prnt = node;
            T.nodes(chldrn(i)).lvl = lvl+1;
            end
        end
    end
    if lvl>0
        for nk=1:4^lvl
            node = nk+n_idx;
            chldrn = T.nodes(node).chld;
            for i=1:length(chldrn)
                T = siblings(T,i,chldrn);
            end
        end
    end
    n_idx = n_idx + 4^lvl;
end

T.n_blocks=size(T.nodes,2);

T.lvp = [T.lvp T.n_blocks];

end

function T=siblings(T,i,chldrn)
  if i==1
     
     T.nodes(chldrn(i)).sbl_r = chldrn(2);
     T.nodes(chldrn(i)).sbl_d = chldrn(3);
     p = T.nodes(chldrn(i)).prnt;
     s1 = T.nodes(p).sbl_u;
     s2 = T.nodes(p).sbl_l;
     if ~isempty(s1)
         T.nodes(chldrn(i)).sbl_u = T.nodes(s1).chld(3);
     else
         T.nodes(chldrn(i)).sbl_u = T.nodes(T.nodes(p).sbl_d).chld(3);
     end
     if ~isempty(s2)
         T.nodes(chldrn(i)).sbl_l = T.nodes(s2).chld(2);
     else
         T.nodes(chldrn(i)).sbl_l = T.nodes(T.nodes(p).sbl_r).chld(2); 
     end
     T.nodes(chldrn(i)).sbls = [chldrn(2), chldrn(3), chldrn(4), ...
         T.nodes(chldrn(i)).sbl_l, T.nodes(chldrn(i)).sbl_u, T.nodes(T.nodes(p).sbl_u).chld(4),...
         T.nodes(T.nodes(p).sbl_l).chld(4),T.nodes(T.nodes(T.nodes(p).sbl_u).sbl_l).chld(4)];
     T.nodes(chldrn(i)).sbls = sort(T.nodes(chldrn(i)).sbls);
  elseif i==2
     T.nodes(chldrn(i)).sbl_l = chldrn(1);
     T.nodes(chldrn(i)).sbl_d = chldrn(4); 
     p = T.nodes(chldrn(i)).prnt;
     s1 = T.nodes(p).sbl_u;
     s2 = T.nodes(p).sbl_r;
     if ~isempty(s1)
         T.nodes(chldrn(i)).sbl_u = T.nodes(s1).chld(4);
     else
         T.nodes(chldrn(i)).sbl_u = T.nodes(T.nodes(p).sbl_d).chld(4);
     end
     if ~isempty(s2)
         T.nodes(chldrn(i)).sbl_r = T.nodes(s2).chld(1);
     else
         T.nodes(chldrn(i)).sbl_r = T.nodes(T.nodes(p).sbl_l).chld(1);
     end
     T.nodes(chldrn(i)).sbls = [T.nodes(chldrn(i)).sbl_d, T.nodes(chldrn(i)).sbl_r, ...
         T.nodes(chldrn(i)).sbl_l, T.nodes(chldrn(i)).sbl_u, T.nodes(T.nodes(p).sbl_u).chld(3),...
         T.nodes(T.nodes(p).sbl_r).chld(3),T.nodes(T.nodes(T.nodes(p).sbl_u).sbl_r).chld(3),...
         chldrn(3)];
     T.nodes(chldrn(i)).sbls = sort(T.nodes(chldrn(i)).sbls);
  elseif i==3
     T.nodes(chldrn(i)).sbl_u = chldrn(1);
     T.nodes(chldrn(i)).sbl_r = chldrn(4);
     p = T.nodes(chldrn(i)).prnt;
     s1 = T.nodes(p).sbl_l;
     s2 = T.nodes(p).sbl_d;
     if ~isempty(s1)
         T.nodes(chldrn(i)).sbl_l = T.nodes(s1).chld(4);
     else
         T.nodes(chldrn(i)).sbl_l = T.nodes(T.nodes(p).sbl_r).chld(4);
     end
     if ~isempty(s2)
         T.nodes(chldrn(i)).sbl_d = T.nodes(s2).chld(1);
     else
         T.nodes(chldrn(i)).sbl_d = T.nodes(T.nodes(p).sbl_u).chld(1);
     end
     T.nodes(chldrn(i)).sbls = [T.nodes(chldrn(i)).sbl_d, T.nodes(chldrn(i)).sbl_r, ...
         T.nodes(chldrn(i)).sbl_l, T.nodes(chldrn(i)).sbl_u, T.nodes(T.nodes(p).sbl_d).chld(2),...
         T.nodes(T.nodes(p).sbl_l).chld(2),T.nodes(T.nodes(T.nodes(p).sbl_d).sbl_l).chld(2),...
         chldrn(2)];
     T.nodes(chldrn(i)).sbls = sort(T.nodes(chldrn(i)).sbls);
  else
     T.nodes(chldrn(i)).sbl_u = chldrn(2);
     T.nodes(chldrn(i)).sbl_l = chldrn(3);
     p = T.nodes(chldrn(i)).prnt;
     s1 = T.nodes(p).sbl_r;
     s2 = T.nodes(p).sbl_d;
     if ~isempty(s1)
         T.nodes(chldrn(i)).sbl_r = T.nodes(s1).chld(3);
     else
         T.nodes(chldrn(i)).sbl_r = T.nodes(T.nodes(p).sbl_l).chld(3);
     end
     if ~isempty(s2)
         T.nodes(chldrn(i)).sbl_d = T.nodes(s2).chld(2);
     else
         T.nodes(chldrn(i)).sbl_d = T.nodes(T.nodes(p).sbl_u).chld(2);
     end
     T.nodes(chldrn(i)).sbls = [T.nodes(chldrn(i)).sbl_d, T.nodes(chldrn(i)).sbl_r, ...
         T.nodes(chldrn(i)).sbl_l, T.nodes(chldrn(i)).sbl_u, T.nodes(T.nodes(p).sbl_d).chld(1),...
         T.nodes(T.nodes(p).sbl_r).chld(1),T.nodes(T.nodes(T.nodes(p).sbl_d).sbl_r).chld(1),...
         chldrn(1)];
     T.nodes(chldrn(i)).sbls = sort(T.nodes(chldrn(i)).sbls);
  end
end