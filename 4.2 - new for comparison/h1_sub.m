% the h1-subproblem and set g_h1.
% h1 = lambda1*||x||_1

function z_h1 = h1_sub(g_f, g_h2, xk, d, lambda1)  

tau = xk - (g_f + g_h2)./diag(d);

tau_sign = tau;
tau_sign(tau>=0) = 1;
tau_sign(tau<0) = -1;
z_h1 = tau_sign.*max(0,abs(tau) - lambda1./diag(d));
%fval = (g_f')*z_h1 + yh1(z_h1, lambda1) + (g_h2')*z_h1 + 1/2*((z_h1 - xk)')*(z_h1-xk);

end



