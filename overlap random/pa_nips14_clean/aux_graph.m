function [grad, ell, f, z] = aux_graph(x, XX, Xy, y, gamma, alpha, mu, E)
grad = gamma * (XX*x - Xy);
if nargout == 1
    return;
end

ell = x'*grad/2 + gamma/2*(y'*y - x'*Xy);

reg = 0;
z = x;
for ii = 1:size(E,1)
    e = E(ii,:);
    dif = abs(diff(x(e)));
    reg = reg + alpha(ii)*dif;
    sgn = sign(x(e(1))-x(e(2)));
    if dif > 2*mu
        z(e) = z(e) - sgn*mu*alpha(ii);
    else
        z(e) = z(e) - sgn*dif/2*alpha(ii);
    end
end
f = ell + reg;