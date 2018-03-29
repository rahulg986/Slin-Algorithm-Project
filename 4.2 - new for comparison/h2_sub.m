% the h2-subproblem and set g_h2.
% h2 = lambda2*|Rx|_1

%% Define the function to solve the subproblem:
function z_h2 = h2_sub(R, g_f, g_h1, xk, d, lambda2)

% This is the main function to compute the h2-subproblem. All the input arguments correspond to equation (33) in our written draft.
p = length(xk);

Q = sparse(R*(d\(R')));
b = -R*(xk - d\(g_f + g_h1));


% Call Minh's function
mu = block_genR(Q, b, zeros(p-1,1)-lambda2, zeros(p-1,1)+lambda2, zeros(p-1,1), 100);

% Figure out the original solution
z_h2 = xk - d\(g_f + g_h1 + (R')*mu);
%fval = (g_f')*z_h2 + (g_h1')*z_h2 + yh2(z_h2, R, lambda2) + 1/2*((z_h2 - xk)')*(z_h2-xk);
end