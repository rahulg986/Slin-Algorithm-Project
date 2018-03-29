% the h2-subproblem and set g_h2.
% h2 = lambda2*|Rx|_1

%% Define the function to solve the subproblem:
function z_h = hi_sub(next_index, a_f, a_h, T, x, d, lambda1)
l = a_f + sum(a_h,2)-a_h(:,next_index);
% This is the main function to compute the h2-subproblem. 
d = diag(d);  %% Make d as a vector!
p = size(a_f,1);
z_h = zeros(p,1);
T_1 = T(next_index,:)==1;
T_0 = T(next_index,:)==0;
%t = sqrt(sumsqr(m));
%k = d(1,1)/(1-lambda1/t);
%z_h(T_1) = m / k;
k_min = 0;
dmax = max(d(T_1));
N = size(d(T_1),1);
c2min = min((d(T_1).*x(T_1) - l(T_1)).^2);
k_max = dmax/(sqrt(N*c2min/lambda1^2)-1);
epsilon = 1e-5;
k_left = k_min;
k_right = k_max;
while (k_right - k_left >= epsilon) 
    k_test = (k_right + k_left)/2;
    m = (d(T_1).*x(T_1) - l(T_1))./(1+d(T_1)/k_test);
    value = sumsqr(m);
if (value <= lambda1^2)
 k_left = k_test;
else  % (value > dg^2)
k_right = k_test;
end
end  % end while
k = (k_right + k_left)/2;
z_h(T_1) = (d(T_1).*x(T_1) - l(T_1))./(k + d(T_1));
z_h(T_0) = x(T_0) - l(T_0)./d(T_0);

% Figure out the original solution
%z_h2 = xk - d\(g_f + g_h1 + (R')*mu);
%fval = (g_f')*z_h2 + (g_h1')*z_h2 + yh2(z_h2, R, lambda2) + 1/2*((z_h2 - xk)')*(z_h2-xk);
end