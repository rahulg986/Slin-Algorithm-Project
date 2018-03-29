% the h1-subproblem and set g_h1.
% h1 = lambda1*||x||_1

function z_h1 = vu_h1_sub(xk, T1, w1)  

tau = T1;

tau_sign = tau;
tau_sign(tau>=0) = 1;
tau_sign(tau<0) = -1;
z_h1 = tau_sign.*max(0,abs(tau) - w1);
%fval = (g_f')*z_h1 + yh1(z_h1, lambda1) + (g_h2')*z_h1 + 1/2*((z_h1 - xk)')*(z_h1-xk);

end



