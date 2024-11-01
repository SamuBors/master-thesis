% OUTPUTS
% r = roots of the determinant equation

% INPUTS
% M = cell (1 x k) of VAR matrices

function [r] = VAR_roots(M)

    syms z

    s = eye(size(M{1,1},1));
    for q = 1:size(M,2)
        s = s-M{1,q}*z^q;
    end
    d = det(s);
    r = double(vpa(root(d,z)));

end