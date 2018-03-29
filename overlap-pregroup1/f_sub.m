%Solve the f-subproblem and set g^k_f.

%% We define the procedure of solving the subproblem as a function:
function z_f = f_sub(A, b, sum_gh, xk, d,L,lamda,INV)

%% This is the main function to compute the f-subproblem. All the input arguments correspond to equation (28) in our written draft.
%method direct inversion
M = d*xk - sum_gh + 1/(lamda*L)*A'*b;
z_f = INV*M;
%method preconjgrad
% p = length(xk);
% 
% x0 = zeros(1,p);
% itmax = 1000;
% t=1.0;
% tol = 1.0e-03;
% V = d;
% input.A = A;
% input.b = b;
% input.sum_gh = sum_gh;
% input.xk = xk;
% input.d = d;
% input.L = L;
% input.lamda = lamda;



%[z_f,~] = PreConjGrad(@ysubf,input,x0,V,t,itmax,tol);
%z_f = z_f';
%N=d;
%z_f = xk+ pcg(1/(L*lamda)*(A')*A + d,(A')*b/(L*lamda) - sum_gh + A'*(A*xk),1e-3,[],N);
%z_f = xk+ pcg(1/(L*lamda)*(A')*A + d,(A')*b/(L*lamda) - sum_gh + d*xk),1e-3,[],N);
%method quadprog

%opts = optimoptions('quadprog', 'Diagnostics', 'off', 'Display', 'off');
%z_f = quadprog(1/(L*lamda)*(A')*A + d, -(A')*b/(L*lamda) + sum_gh - d*xk, [], [], [], [], [], [], [], opts);

end


