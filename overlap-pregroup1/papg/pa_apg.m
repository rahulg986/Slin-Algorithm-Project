function [y, f] = pa_apg(x, mu, aux, maxiter)

y = x;
f = zeros(maxiter, 1);
eta = 1;
for ii = 1:maxiter
    x_old = x;
    %[grad,~,f(ii),~] = aux(y);
    z = y - mu*aux(y);
    [~,~,f(ii),x] = aux(z);
    eta_old = eta;
    eta = (1+sqrt(1+4*eta^2)) / 2;
    y = x + (eta_old-1)/eta*(x-x_old);
    if mod(ii,100) == 0
        fprintf('pa_apg: iter = %d, f = %f \n', ii, f(ii));
    end
end
