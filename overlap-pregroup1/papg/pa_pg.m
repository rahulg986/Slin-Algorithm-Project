function [x, f] = pa_pg(x, mu, aux, maxiter)

f = zeros(maxiter, 1);
for ii = 1:maxiter
    %[grad,~,f(ii),~] = aux(x);
    z = x - mu*aux(x);
    [~,~,f(ii),x] = aux(z);
    if mod(ii,100) == 0
        fprintf('pa_pg: iter = %d, f = %f \n', ii, f(ii));
    end    
end