function history = pa_pg(x, mu, aux, maxiter)

t_start = tic;
f = zeros(maxiter, 1);
for ii = 1:maxiter
    %[grad,~,f(ii),~] = aux(x);
    z = x - mu*aux(x);
    xold = x;
    [~,~,f(ii),x] = aux(z);
    
    history.time(ii) = toc(t_start);
    history.objval(ii) = f(ii);
    history.res(ii) = norm(x-xold)/norm(xold);
    
    %fprintf('pa_pg: iter = %d, f = %f, res = %f, time = %f \n', ii, f(ii), history.res(ii), history.time(ii));
    
    if history.res(ii) < 1e-4
        break;
    end
    
%     if mod(ii,100) == 0
%         fprintf('pa_pg: iter = %d, f = %f \n', ii, f(ii));
%     end    
end