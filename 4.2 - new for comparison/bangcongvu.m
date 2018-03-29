%primal-dual Bang cong vu
%% The objective function is: min 1/2 |Ax-b|^2 + lambda_1 |x|_1 + lambda_2 sum_{j} |x_{j+1} - x_{j}|_1
function [time,x,iter_num,optimal_value_mark,optimal_value] = bangcongvu(input)
% Same use of time_limit here!

warning('off', 'all');

%% Get the data from "input".
n = input.n;
p = input.p;
A = input.A;
b = input.b;
% Precompute the R function
R = zeros(p-1,p);
for i = 1:(p-1)
    R(i,i) = -1;
    R(i,i+1) = 1;
end
L2=R;
L1 = eye(p);
w1 = input.lambda1; %/(input.lambda1 + input.lambda2);
w2 = input.lambda2; %/(input.lambda1 + input.lambda2);
tao=10/norm(diag(A'*A));
%tao = 0.05;
EPI = input.EPI; 

%% Data preparation
MAX_ITER = 1000;

% Initial prox center
x = ones(p,1);  % By default, we assume everything to be column vector.

% Inital values of h1 and h2 subproblems.
S= x;
T1=x;
T2=ones(p-1,1);
v1=x;
v2=ones(p-1,1);
pn=x;
z_h2=ones(p-1,1);
%z_h1 = x;
%z_h2 = x;

%g_h1 = gh1(z_h1, lambda1);
%g_h2 = gh2(z_h2, R, lambda2);

%beta1 = 0.5;

% Prepare the output structure
optimal_value = zeros(MAX_ITER,1);  % We only record the cases where the prox center changes!
optimal_value(1) = objfunc(x, A, b, w1, L2, w2);
optimal_value_mark = 1;


%fprintf('%10s\t%10s\n', optimal_value(1),  v_values(1));

iter_num = 0;

%%%algorithm start
tic;   %start time
for iter_num = 1:MAX_ITER
    % Check whether the stuck is from ALIN or from quadratic minimization!
    % iter_num = iter_num + 1;      
    if (iter_num > MAX_ITER)
        save('bad_matrix.mat', 'A', 'b');
        break;        
    end

    % Solve f-subproblem and set the corresponding g.
        S= x - tao*(w1*L1'*v1 + w2*L2'*v2);
        pn = vu_f_sub(A, b, x, S);
        yn=2*pn-x;
        x = pn;
        optimal_value_mark = optimal_value_mark + 1;
        optimal_value(optimal_value_mark) = objfunc(x, A, b, w1, L2, w2);
        %solve h1-sub
        T1= v1 + L1*yn;
        q1n = v1- vu_h1_sub(x, T1,w1);
        v1=q1n; 
        %solve h2-sub
        T2= v2 + L2*yn;
        q2n = v2- vu_h2_sub(L2(1:p-2,1:p-1), z_h2, T2, w2);
        v2=q2n;
        fprintf('optimal value = %10s\n', optimal_value(optimal_value_mark));
        time = toc;
%         if (time > time_limit)
%             break;
%         end
end

% Output the information for this iteration!
%fprintf('The %d iteration\t', iter_num);
%fprintf('%10s\n', optimal_value(optimal_value_mark));

time = toc;
%time = toc;
fprintf('The total number of iterations is %d\n', iter_num);
%fprintf('The objective value is %d\n', optimal_value);
% figure
% %Plot the optimal value changes with respect to the number of iterations
% plot(1:optimal_value_mark,optimal_value(1:optimal_value_mark),'-*');
% ylabel('optimal value');
% xlabel('Iterations');
% title('Change of optimal value w.r.t. the number of iterations');

% Plot the change of ln(F - F^*) with respect to the number of iterations
%  figure
%  value_diff = optimal_value(1:optimal_value_mark) - optimal_value(optimal_value_mark);
% % if (update_selector ==  1)
% %     plot(1:(optimal_value_mark-1), log(value_diff(1:optimal_value_mark-1)), '-*');
% %     hold on;
% % else
% %     plot(1:(optimal_value_mark-1), log(value_diff(1:optimal_value_mark-1)), '-o');
% % end
%  plot(1:(optimal_value_mark-1), log(value_diff(1:optimal_value_mark-1)), '-*');
% xlabel('Iterations');
% ylabel('ln(F(x) - F(x*))');
% title('Change of ln(F(x) - F(x^*)) w.r.t. the number of iterations');
end
