function history = s_pg(x, mu, L0, aux, maxiter)

t_start = tic;
%f = zeros(maxiter, 1);
beta = 1 / (1+mu*L0);
for ii = 1:maxiter
    [grad,~,f(ii),y] = aux(x);
    z = x - grad/L0;
    xold = x;
    x = (1-beta)*z + beta*y;
    
    history.time(ii) = toc(t_start);
    history.objval(ii) = f(ii);
    history.res(ii) = norm(x-xold);
    
    %fprintf('s_pg: iter = %d, f = %f, res = %f, time = %f \n', ii, f(ii), history.res(ii), history.time(ii));
    
    if history.res(ii) < 1e-4
        break;
    end
%     if mod(ii,100) == 0
%         fprintf('s_pg: iter = %d, f = %f \n', ii, f(ii));
%     end    
end