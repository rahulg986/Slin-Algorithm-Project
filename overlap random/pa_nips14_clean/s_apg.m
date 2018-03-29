function history = s_apg(x, mu, L0, aux, maxiter)

t_start = tic;
y = x;
f = zeros(maxiter, 1);
eta = 1;
%y = 0;
for ii = 1:maxiter
    x_old = x;
    [grad,~,f(ii),z] = aux(y);
    beta = 1 / (1+mu*L0);
    x = (1-beta)*(y - grad/L0) + beta*z;
    eta_old = eta;
    eta = (1+sqrt(1+4*eta^2)) / 2;
    yold = y;
    y = x + (eta_old-1)/eta*(x-x_old);
    
    history.time(ii) = toc(t_start);
    history.objval(ii) = f(ii);
    history.res(ii) = norm(y-yold)/norm(yold);
    
    %fprintf('pa_pg: iter = %d, f = %f, res = %f, time = %f \n', ii, f(ii), history.res(ii), history.time(ii));
    
    if history.res(ii) < 1e-4
        break;
    end
    
%     if mod(ii,100) == 0
%         fprintf('s_apg: iter = %d, f = %f \n', ii, f(ii));
%     end
end