function [x, f] = s_pg(x, mu, L0, aux, maxiter)

f = zeros(maxiter, 1);
beta = 1 / (1+mu*L0);
for ii = 1:maxiter
    [grad,~,f(ii),y] = aux(x);
    z = x - grad/L0;
    x = (1-beta)*z + beta*y;
    if mod(ii,100) == 0
        fprintf('s_pg: iter = %d, f = %f \n', ii, f(ii));
    end    
end