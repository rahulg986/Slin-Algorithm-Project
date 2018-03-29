%Solve the f-subproblem and set g^k_f.

%% We define the procedure of solving the subproblem as a function:
function z_f = f_sub(A, b, sum_gh, xk, d,L,lamda,INV, z)
%tic
%% This is the main function to compute the f-subproblem. All the input arguments correspond to equation (28) in our written draft.
%method direct inversion
M = d*xk - sum_gh + 1/(lamda*L)*(A')*b;

%% Use Minh's code
[z_f, ~] = prec_conjug(1/sqrt(lamda*L)*A, M, diag(d), z);
%z_f = INV\M;

%z_f = (LA') \ (LA \ M);
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

%matlab builtin PCG function

%z_f = pcg(INV,M,[],[],inv(d));

%method quadprog

%opts = optimoptions('quadprog', 'Diagnostics', 'off', 'Display', 'off');
%z_f = quadprog(1/(L*lamda)*(A')*A + d, -(A')*b/(L*lamda) + sum_gh - d*xk, [], [], [], [], [], [], [], opts);


%toc
end


