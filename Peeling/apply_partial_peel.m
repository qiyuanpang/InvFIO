function [F] = apply_partial_peel(GHat,T,Y,x,L)
N = size(x,2);

m = size(Y,2);


 % Initialize with zeros before we tack things on
 F = zeros(N,m);
%  iLevel_threshold = max(T.nlvl,4);
iLevel_threshold = 4;
 for iPrev = iLevel_threshold:L
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
     
end
