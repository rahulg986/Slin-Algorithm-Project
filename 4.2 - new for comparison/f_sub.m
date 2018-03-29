%Solve the f-subproblem and set g^k_f.

%% We define the procedure of solving the subproblem as a function:
function z_f = f_sub(A, b, g_h1, g_h2, xk, d)

%% This is the main function to compute the f-subproblem. All the input arguments correspond to equation (28) in our written draft.

p = length(xk);

x0 = zeros(1,p);
itmax = 20;
t=1.0;
tol = 1.0e-03;
V = diag(d);

%using quadprog method
%opts = optimoptions('quadprog', 'Diagnostics', 'off', 'Display', 'off');
%[z_f,~] = PreConjGrad(@ysubf,input,x0,V,t,itmax,tol);

%z_f = quadprog((A')*A + d, -(A')*b + g_h1 + g_h2 - d*xk, [], [], [], [], [], [], [], opts);
%using pcg
M = d;
%z_f = pcg((A')*A + d, (A')*b - g_h1 - g_h2 + d*xk);
%z_f = pcg((A')*A + d, (A')*b - g_h1 - g_h2 + d*xk,[],[],M);
z_f = xk+ pcg(((A')*A + d),(A')*b - g_h1 - g_h2 - A'*(A*xk),1e-3,[],M);
% %using minh's conjugate
% x_cur = xk;
% c = diag(A'*A);
% z_f = prec_conjug(A,(A')*b - g_h1 - g_h2 + d*xk,c,x_cur);

%fval = 1/2*(norm(A*z_f-b)^2) + ((g_h1+g_h2)')*z_f + 1/2*((z_f - xk)')*(z_f-xk);

end


