function [y, f] = s_apg(x, mu, L0, aux, maxiter)

y = x;
f = zeros(maxiter, 1);
eta = 1;
for ii = 1:maxiter
    x_old = x;
    [grad,~,f(ii),z] = aux(y);
    beta = 1 / (1+mu*L0);
    x = (1-beta)*(y - grad/L0) + beta*z;
    eta_old = eta;
    eta = (1+sqrt(1+4*eta^2)) / 2;
    y = x + (eta_old-1)/eta*(x-x_old);
    if mod(ii,100) == 0
        fprintf('s_apg: iter = %d, f = %f \n', ii, f(ii));
    end
end