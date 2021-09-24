%    References:
%
%      Lin Lin, Jianfeng Lu, and Lexing Ying. 2011. Fast construction of
%      hierarchical matrix representation from matrix-vector multiplication.
%      J. Comput. Phys. 230, 10 (May 2011), 4071-4087.

function [F] = apply_peel(GHat,T,Y,x)

N = size(x,2);

m = size(Y,2);

nLevel = T.nlvl;

% Initialize with zeros before we tack things on
F = zeros(N,m);
% iLevel_threshold = max(T.nlvl,4);
iLevel_threshold = 4;
for iPrev = iLevel_threshold:nLevel
    n_nodes = T.lvp(iPrev+1) - T.lvp(iPrev);

    % Compute transfer
    for I = 1:n_nodes
        iNodeIdx = T.lvp(iPrev) + I;
        inter = T.nodes(iNodeIdx).inter;
        inter = inter - T.lvp(iPrev);
        for J = inter
            if I < J
                M_IJ = GHat{iPrev}(I,J).M;
            else
                M_IJ = GHat{iPrev}(J,I).M';
            end
            jNodeIdx = T.lvp(iPrev) + J;
            xj       = T.nodes(jNodeIdx).xi;
            xi = T.nodes(iNodeIdx).xi;
            
            F(xi,:) = F(xi,:) + GHat{iPrev}(I,J).U*(M_IJ * (GHat{iPrev}(J,I).U'*Y(xj,:)));
        end
    end
end

n_nodes = T.lvp(nLevel+1) - T.lvp(nLevel);
for I = 1:n_nodes
    iNodeIdx = T.lvp(nLevel) + I;
    xi = T.nodes(iNodeIdx).xi;
    F(xi,:) = F(xi,:) + GHat{nLevel+1}(I,I).M*Y(xi,:);
    for jNodeIdx = T.nodes(iNodeIdx).nbor 
        if jNodeIdx~=iNodeIdx
            xj = T.nodes(jNodeIdx).xi;
            J = jNodeIdx - T.lvp(nLevel);
            F(xi,:) = F(xi,:) + GHat{nLevel+1}(I,J).M*Y(xj,:);
        end
    end
end

end
