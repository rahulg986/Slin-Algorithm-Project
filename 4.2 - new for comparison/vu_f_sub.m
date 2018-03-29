%Solve the f-subproblem and set g^k_f.

%% We define the procedure of solving the subproblem as a function:
function z_f = vu_f_sub(A, b, xk, s)

%% This is the main function to compute the f-subproblem. All the input arguments correspond to equation (28) in our written draft.

p = length(xk);

x0 = zeros(1,p);
itmax = 1000;
t=1.0;
tol = 1.0e-03;
%V = diag(d);
 
opts = optimoptions('quadprog', 'Diagnostics', 'off', 'Display', 'off');
z_f = quadprog((A')*A + eye(p), -(A')*b + s, [], [], [], [], [], [], [],opts);%my understanding
%z_f = quadprog((A')*A + norm(diag(A'*A))*eye(p), -(A')*b +
%norm(diag(A'*A))*s, [], [], [], [], [], [], [], opts);%R ge Understanding
%-> not converge
%fval = 1/2*(norm(A*z_f-b)^2) + ((g_h1+g_h2)')*z_f + 1/2*((z_f - xk)')*(z_f-xk);

end


