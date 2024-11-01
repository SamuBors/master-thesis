% OUTPUTS
% r = roots of the determinant equation

% INPUTS
% M = cell (1 x q) of MA matrices

function [r] = VMA_roots(M)
    
    if size(M,2)>1
        syms z
        
        s = M{1,1};
        for q = 2:size(M,2)
            s = s+M{1,q}*z^(q-1);
        end
        d = det(s);
        r = double(vpa(root(d,z)));
    else
        r = 'Not a VMA';
    end

end