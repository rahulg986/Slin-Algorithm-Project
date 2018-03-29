function [grad, ell, f, z] = aux_group(x, X, Xy, y, gamma, alpha, mu, groups)
Xy = X'*y;
grad = gamma * ( X'*(X*x) - Xy);
if nargout == 1
    return;
end

ell = x'*grad/2 + gamma/2*(y'*y - x'*Xy);

reg = 0;
z = x;
for ii = 1:length(groups)
    group = groups{ii};
    normg = norm(x(group),2);
    reg = reg + alpha(ii)*normg;
    
    if normg > mu
        z(group) = z(group) - mu/normg*x(group)*alpha(ii);
    else
        z(group) = z(group) - x(group)*alpha(ii);
    end
end
f = ell + reg;