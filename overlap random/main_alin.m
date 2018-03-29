%% The objective function is: min 1/(2*L*lamda) |Ax-b|^2 + sum_g d_g ||x_g||_2
function [time,x,iter_num,optimal_value_mark,optimal_value] = main_alin(input, update_selector,p_opt)

warning('off', 'all');
%% Get the data from "input".
n = input.n;
p = input.p;
A = input.A;
b = input.b;
group_ind=input.group_ind;
lambda1 = input.lambda1;%1/100
% mean_d = mean(diag(A'*A));
% d = zeros(p,1) + mean_d;
% d = diag(d);
%d = 0.0001;

%EPI = input.EPI; 
EPI = 1e-4;
T = input.T;
L = input.ng;
lamda = L/5;
ng = input.ng;
d = 1/(ng*lamda)*diag(diag(A'*A));
%d = 10/(ng*lamda)*diag(diag(A'*A));
%d = diag(diag(A'*A));
%% Data preparation
MAX_ITER = 2000;
%% Comput the inverse of matrix
%INV = inv(1/(lamda*L)*(A'*A) + d(1,1)*eye(p));  %% Note: this is for the
%case of a single rho.
INV = 1/(lamda*L)*(A'*A) + d;  %% This is for general diagonal matrix d.
% LU Decomposition
%LA = chol(INV, 'lower');
%LA = sparse(LA);
% Initial prox center
x = ones(p,1);  % By default, we assume everything to be column vector.
z_f = zeros(p,1);
a_f = zeros(p,1);
b_f = 0;
z_h = zeros(p,ng);
a_h = zeros(p,ng);
b_h = zeros(1,ng);
% Inital values of hi subproblems.
for i = 1:ng
z_h(:,i) = x;
a_h(:,i) = gh(z_h(:,i), lambda1,T(i,:));
b_h(i)= yhi(i,z_h(:,i),lambda1,group_ind) - a_h(:,i)'*z_h(:,i);
end

beta1 = 0.2;

% Prepare the output structure
optimal_value = zeros(MAX_ITER,1);  % We only record the cases where the prox center changes!
optimal_value(1) = objfunc(x, A, b, lamda, L,lambda1,T);
optimal_value_mark = 1;

% v values 
v_values = zeros(MAX_ITER,1);  % Record the v values
v_values(1) = 1;
v_values_mark = 1;

fprintf('%10s\t%d\n', optimal_value(1),  v_values(1));

iter_num = 0;

%% update_selector: 0: No update rule nor selection rule; 1: Has update rule but without selection rule; 2: Both rules are included.
tic;  %start time
next_index = 0;
%while (v_values(v_values_mark) >= EPI)
while (true)
    % Check whether the stuck is from SLIN or from quadratic minimization!
    iter_num = iter_num + 1;      
    if (iter_num > MAX_ITER)
        save('bad_matrix.mat', 'A', 'b');
        break;        
    end

    % Solve the (next_index)-subproblem and set the corresponding g.
    if (next_index == 0)
        sum_ah = sum(a_h,2);
        %% Solving the f-subproblem with our own code        
        z_f = f_sub(A, b, sum_ah, x, d, L, lamda,INV, z_f);
        %% Verify the computed value is really optimal
        if (max(abs(INV*z_f - (d*x - sum_ah + 1/(lamda*L)*(A')*b))) > 1e-5)
            fprintf('Stop here!');
        end
        %% Solving the f-subproblem with MATLAB fminunc
        %fsub_fun = @(x)ysubf(x, A, b, sum_ah, x, d, L, lamda);
        %z_f0 = ones(p,1);
        %[z_f, ~] = fminunc(fsub_fun, z_f0);
        a_f = -sum_ah - d*(z_f - x);
        b_f = yf(z_f,A, b,L,lamda) - a_f'*z_f;     
        temp_x = z_f;
        temp_ftilde = yf(temp_x, A, b, L, lamda) + sum(a_h'*temp_x) + sum(b_h);
        
    else
        %% Solving the hi-subproblem with our own code
        [z_h(:,next_index), new_d] = hi_sub(next_index, a_f, a_h, T, x, d, lambda1);
        %% Solving the hi-subproblem with MATLAB fminunc
        %hisub_fun = @(x)ysubhi(x, next_index, a_f, a_h, T, x, d, lambda1);
        %z_h0 = ones(p,1);
        %[z_h(:, next_index), ~] = fminunc(hisub_fun, z_h0);
        a_h(:,next_index) = - a_f - (sum(a_h,2)-a_h(:,next_index)) - diag(new_d)*(z_h(:,next_index) - x);
        b_h(next_index) =yhi(next_index,z_h(:,next_index),lambda1,group_ind) - a_h(:,next_index)'*z_h(:,next_index);
        temp_x = z_h(:,next_index);
        temp_ftilde = yhi(next_index,z_h(:,next_index),lambda1,group_ind) + a_f'*temp_x + b_f + (sum(a_h,2)-a_h(:,next_index))'*temp_x + sum(b_h) - b_h(next_index);      
     end

    % Update x^{k+1} and find the next subproblem to solve       
    v = optimal_value(optimal_value_mark) - temp_ftilde;
    v_values_mark = v_values_mark + 1;
    v_values(v_values_mark) = v;
    if (v < 0)
        fprintf('v is negative!');
    end

    % Determine the index of the next sub-problem to solve
    if (update_selector ==  0)  %opterator splitting method
        next_index = mod(next_index + 1, ng+1);
        x = temp_x;
        optimal_value_mark = optimal_value_mark + 1;
        optimal_value(optimal_value_mark) = objfunc(x,  A, b, lamda, L,lambda1,T);
    else    
        curr_index = next_index;
        
        [next_index, total_diff] = max_difference(temp_x, a_f, b_f, a_h, b_h,T, group_ind, A, b,lamda, L,lambda1, next_index);
        
        if (update_selector == 1)
            % Remove the selection rule
            next_index = mod(curr_index + 1, ng+1);           
        end

        if (total_diff <= (1-beta1)*v)            
            % Update the prox center vector and the corresponding optimal
            % value
            x = temp_x;
            optimal_value_mark = optimal_value_mark + 1;
            optimal_value(optimal_value_mark) = objfunc(x, A, b, lamda, L, lambda1, T);                 
%             %for comparison with method
            if (optimal_value(optimal_value_mark)-p_opt <1e-4)
                break;
            end
        else
            % Keep the original prox center
            optimal_value_mark = optimal_value_mark + 1;
            optimal_value(optimal_value_mark) = optimal_value(optimal_value_mark - 1);
        end      
    end

    % Output the information for this iteration!
    fprintf('The %d iteration\t', iter_num);
    fprintf('%10s\t%10s\t next iteration is: %d \n', optimal_value(optimal_value_mark),  v_values(v_values_mark), next_index);
end
toc;
time = toc;
fprintf('The total number of iterations is %d\n', iter_num);
%fprintf('The objective value is %d\n', optimal_value);
% figure
% %Plot the optimal value changes with respect to the number of iterations
% plot(1:optimal_value_mark,optimal_value(1:optimal_value_mark),'-*');
% ylabel('objective value');
% xlabel('Iterations');
% title('Change of optimal value w.r.t. the number of iterations');

% Plot the change of ln(F - F^*) with respect to the number of iterations
% figure
% value_diff = optimal_value(1:optimal_value_mark) - optimal_value(optimal_value_mark);
% if (update_selector ==  1)
%     plot(1:(optimal_value_mark-1), log(value_diff(1:optimal_value_mark-1)), '-*');
%     hold on;
% else
%     plot(1:(optimal_value_mark-1), log(value_diff(1:optimal_value_mark-1)), '-o');
% end
% %plot(1:(optimal_value_mark-1), log(value_diff(1:optimal_value_mark-1)), '-*');
% xlabel('Iterations');
% ylabel('ln(F(x) - F(x*))');
% title('Change of ln(F(x) - F(x^*)) w.r.t. the number of iterations');
end
