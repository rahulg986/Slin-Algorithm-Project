%% The objective function is: min 1/2 |Ax-b|^2 + lambda_1 |x|_1 + lambda_2 sum_{j} |x_{j+1} - x_{j}|_1
function [time, v_values, v_values_mark, x, iter_num, optimal_value_mark, optimal_value] = main_alin(input, update_selector)
% Cheng added: add the time_limit argument, so that if the algorithm's run
% time is longer than the time_limit, we force the algorithm to quit. This
% is to deal with the issue of non-convergence!

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
lambda1 = input.lambda1;
lambda2 = input.lambda2;
if(update_selector == -4)
    mean_d=2*mean(mean(diag(diag(A'*A))));
    d = zeros(p,1) + mean_d;
    d = diag(d);
else d=diag(diag(A'*A));
end
EPI = input.EPI; 

%% Data preparation
MAX_ITER = 1000;

% Initial prox center
x = ones(p,1);  % By default, we assume everything to be column vector.

% Inital values of h1 and h2 subproblems.
z_h1 = x;
z_h2 = x;

g_h1 = gh1(z_h1, lambda1);
g_h2 = gh2(z_h2, R, lambda2);

beta1 = 0.5;

% Prepare the output structure
optimal_value = zeros(MAX_ITER,1);  % We only record the cases where the prox center changes!
optimal_value(1) = objfunc(x, A, b, lambda1, R, lambda2);
optimal_value_mark = 1;

% v values 
v_values = zeros(MAX_ITER,1);  % Record the v values
v_values(1) = 1;
v_values_mark = 1;

fprintf('%10s\t%10s\n', optimal_value(1),  v_values(1));

iter_num = 0;

%% update_selector: 0: No update rule nor selection rule; 1: Has update rule but without selection rule; 2: Both rules are included.
tic;   %start time
next_index = 0;
while (v_values(v_values_mark) >= EPI)
    % Check whether the stuck is from ALIN or from quadratic minimization!
    iter_num = iter_num + 1;      
    if (iter_num > MAX_ITER)
        save('bad_matrix.mat', 'A', 'b');
        break;        
    end

    % Solve the (next_index)-subproblem and set the corresponding g.
    if (next_index == 0)
        z_f = f_sub(A, b, g_h1, g_h2, x, d);
        g_f = -g_h1 - g_h2 - d*(z_f - x);
        temp_x = z_f;
        temp_ftilde = yf(temp_x, A, b) + yh1(z_h1, lambda1) + (g_h1')*(temp_x - z_h1) + yh2(z_h2, R, lambda2) + (g_h2')*(temp_x - z_h2);
    elseif (next_index == 1)
        z_h1 = h1_sub(g_f, g_h2, x, d, lambda1);
        g_h1 = -g_f - g_h2 - d*(z_h1 - x);
        temp_x = z_h1;
        temp_ftilde = yf(z_f, A, b) + (g_f')*(temp_x - z_f) + yh1(temp_x, lambda1) + yh2(z_h2, R, lambda2) + (g_h2')*(temp_x - z_h2);
    else
        z_h2 = h2_sub(R, g_f, g_h1, x, d, lambda2);
        g_h2 = -g_f - g_h1 - d*(z_h2 - x);
        temp_x = z_h2;
        temp_ftilde = yf(z_f, A, b) + (g_f')*(temp_x - z_f) + yh1(z_h1, lambda1) + (g_h1')*(temp_x - z_h1) + yh2(temp_x, R, lambda2);
    end

    % Update x^{k+1} and find the next subproblem to solve       
    v = optimal_value(optimal_value_mark) - temp_ftilde;
    v_values_mark = v_values_mark + 1;
    v_values(v_values_mark) = v;

    % Determine the index of the next sub-problem to solve
    if (update_selector ==  0)
        next_index = mod(next_index + 1, 3);
        x = temp_x;
        optimal_value_mark = optimal_value_mark + 1;
        optimal_value(optimal_value_mark) = objfunc(x, A, b, lambda1, R, lambda2);
    else    
        curr_index = next_index;
        [next_index, total_diff] = max_difference(temp_x, z_f, z_h1, z_h2, g_f, g_h1, g_h2, A, b, lambda1, R, lambda2, next_index, 1);
        if (update_selector == 1 ||update_selector == -1||update_selector ==-2||update_selector==-3||update_selector==-4)
            % Remove the selection rule
            next_index = mod(curr_index + 1, 3); 
        end
    if(update_selector > -1) %%for methods slin and alin
        if (total_diff <= (1-beta1)*v)            
            % Update the prox center vector and the corresponding optimal
            % value
            x = temp_x;
            optimal_value_mark = optimal_value_mark + 1;
            optimal_value(optimal_value_mark) = objfunc(x, A, b, lambda1, R, lambda2);                 
        else
            % Keep the original prox center
            optimal_value_mark = optimal_value_mark + 1;
            optimal_value(optimal_value_mark) = optimal_value(optimal_value_mark - 1);
        end  
    elseif(next_index ==0)
        if(update_selector == -2) %%for douglas over relax form
            x = -0.5*x + 1.5*temp_x;
        elseif(update_selector==-3)
            x = 0.5*x + 0.5*temp_x; %for douglas under relax
        else  x = temp_x; %for douglas normal        
        end
            
            
            %%for douglas reachford method.
    end
    end
   if((update_selector == -2)||(update_selector == -1)||(update_selector==-3||(update_selector==-4)))
    optimal_value_mark = optimal_value_mark + 1;
            optimal_value(optimal_value_mark) = objfunc(x, A, b, lambda1, R, lambda2);   
   end
%     time = toc;
%     if (time_limit >= 0)
%         if (time > time_limit)
%             break;
%         end
%    end
    
    % Output the information for this iteration!
    fprintf('The %d iteration\t', iter_num);
    fprintf('%10s\t%10s\n', optimal_value(optimal_value_mark),  v_values(v_values_mark));
end
toc;
time = toc;
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
