%comparison slin with spg?PDMM, SADMM, PA-APG
clear;
clc;

%% Prepare the data
ng=100;  % number of groups 
g_size=100;  % group size
overlap=10;  % number of variables overlapped between two consecutive groups
n=1000;      % sample size
p=ng*(g_size-overlap)+overlap;   % total number of variables 
input.EPI = 1e-4;
input.n = n;
input.p = p;
input.ng = ng;
input.g_size = g_size;
input.overlap = overlap;

%generate toy data
% A: design matrix
% b: output
% T: ng \times p, indicate each group contains which variables 
% Tx: weight for each group
% x: true regression coefficients

 r_T=kron(1:ng, ones(1,g_size));
    c_T=zeros(ng*g_size, 1);
    s=1;
    for g=1:ng
        c_T((g-1)*g_size+1:g*g_size,1)=s:s+g_size-1;
        s=s+g_size-overlap;
    end
    T=sparse(r_T, c_T, 1, ng, p);
     Tx=ones(ng,1);  % uniform weight
     tmp=1:p;
    x=((-1).^tmp(:)).*(exp(-(tmp(:)-1)/100));  
    sn2 = 1;  % signal to noise ratio
    
    A =randn(n,p);
    b= A*x+sn2*rand(n,1);     
input.A=A;
input.b=b;

input.T = T;
%tau=0.1*norm(A'*b,Inf);
lamda1 = 1/g_size;
input.lambda1 = lamda1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

option.maxiter=10000;  
% max iteration 

option.mu=1e-3;   
% Smoothing Parameter
% Carefully choose mu, a larger mu may not converge while a smaller one 
% leads to slow convergence

option.verbose=true;    
option.display_iter=10;   
option.tol=1e-8;          % tolerance

[X, Y, T, Tw, w]=gentoy_group(n, ng, g_size, overlap); 
% generate toy data
% X: design matrix
% Y: output
% T: 0/1 sparse matrix: # of groups by # of features,
% indicate each group contains which variables 
% Tw: # of groups by 1 column vector,  weight for each group
% w: true regression coefficients

[C, g_idx, CNorm] = pre_group(T, Tw); 
% C: C matrix 
% g_idx: group index for the rows of C
% CNorm: 


 
alpha = 1/ng * ones(ng,1);
gamma1 = ng / 5; % /20
gamma=ng/10;   % regularization parameter for group penalty 
lambda=ng/10;  % regularization parameter for L1-norm
AA = A*A';
L0 = abs(eigs(AA,1))/gamma/ng;
eps = 2/2/L0;
mu = min(1/L0, 2*eps);
Ab = A'*b;              % A'*b
[groups] = gen_groups(ng, g_size, overlap);
x = randn(p, 1);
maxiter = 2000;
prob='group';
 aux_pa = @(x) aux_group(x, A, Ab, b, 1/gamma/ng, alpha, mu, groups);
 aux_s = @(x) aux_group(x, A, Ab, b, 1/gamma/ng, alpha, 2*eps, groups);
% [grad_beta,grad_obj,grad_density,grad_iter,grad_time] = ...
%              SPG(prob, Y, X, gamma, lambda, C, CNorm, option, g_idx);  
% SPG with a pre-compuated Lipschitz constant for medium-scale problems
%figure
%[grad_beta,grad_obj,grad_density,grad_iter,grad_time, L] = ...
%              SPG_linesearch(prob, Y, X, gamma, lambda, C,option, g_idx);  
% SPG with line search on Lipschitz constant for large-scale problems

%% Call the three methods for comparison
% [x0,iter_num0,optimal_value_mark0,optimal_value0] = main_alin(input, 0); 
% hold onb_h(2)= b_h(2) + yhi(i,z_h(:,2),lambda1,T) - a_h(:,2)'*z_h(:,2);
%[x1,iter_num1,optimal_value_mark1,optimal_value1] = main_alin(input, 1);
%  hold on
% value_diff1 = optimal_value1(1:optimal_value_mark1) - optimal_value1(optimal_value_mark1);
[x2,iter_num2,optimal_value_mark2,optimal_value2] = main_alin(input, 2);
history = rpdmm_overlap_for(A,b,ng+1,.1,ng,groups,g_size,overlap,alpha,gamma1);
% %% Averaging to get the strict feasible solution
% hx = history.x;
% hz = history.z;
% ave_x = zeros(size(hz));
% % Group 1
% ave_x(1:90) = (hz(1:90) + hx(1:90))/2;
% ave_x(91:100) = hz(91:100) + hx(91:100);
% % Group 2 to 99
% for k = 2:99
%     b_k = 90*k + 10 - 100 + 1;
%     e_k = 90*k + 10;
%     x_b_k = 100*(k-1)+1;
%     x_e_k = 100*k;
%     ave_x(b_k : b_k + 9) = (ave_x(b_k: b_k + 9) + hx(x_b_k: x_b_k + 9))/3;
%     ave_x(b_k + 10: e_k - 10) = (hz(b_k + 10:e_k - 10) + hx(x_b_k + 10:x_e_k - 10))/2;
%     ave_x(e_k - 9: e_k) = hz(e_k-9: e_k) + hx(x_e_k-9:x_e_k);
% end
% % Group 100
% k = 100;
% b_k = 90*k + 10 - 100 + 1;
% e_k = 90*k + 10;
% x_b_k = 100*(k-1)+1;
% x_e_k = 100*k;
% ave_x(b_k : b_k + 9) = (ave_x(b_k: b_k + 9) + hx(x_b_k: x_b_k + 9))/3;
% ave_x(b_k + 10:e_k) = (hz(b_k+10:e_k) + hx(x_b_k + 10:x_e_k))/2;
% % Figure out the objective value of ave_x
% ave_obj = objfunc(ave_x, input.A, input.b, input.g_size/5, input.g_size, input.lambda1, input.T);

%history1 = admm_overlap_for(A,b,ng+1,.1,ng,groups,g_size,overlap,alpha,gamma1);
 history3 = pa_apg(x, mu, aux_pa, maxiter);
% 
 history4 = s_apg(x, 2*eps, L0, aux_s, maxiter);
 history2 = sadmm_overlap(A,b,ng+1,.1,ng,groups,alpha,gamma1);
%value_diff2 = optimal_value2(1:optimal_value_mark2) - optimal_value2(optimal_value_mark2);

%Plot the optimal value changes with respect to the number of iterations

% figure for blow up
% plot(1:optimal_value_mark0,optimal_value0(1:optimal_value_mark0),'-^');
% hold on;
% %plot(1:optimal_value_mark1,optimal_value1(1:optimal_value_mark1),'-*');
% %hold on;
loglog(1:optimal_value_mark2,optimal_value2(1:optimal_value_mark2));
ylabel('Optimal value');
xlabel('Iterations');
hold on
%loglog(1:grad_iter, grad_obj(1:grad_iter), '-^');
%hold on
loglog(1:length(history.objval), history.objval);
hold on
%loglog(1:length(history1.objval), history1.objval, '-+');
%hold on
loglog(1:length(history3.objval), history3.objval);
hold on
loglog(1:length(history4.objval), history4.objval);
hold on
loglog(1:length(history2.objval), history2.objval);
% title('Change of optimal value w.r.t. the number of iterations');
%%%%%%%%figure ln(f(x)-f(x*))
% figure
% %value_diff = optimal_value(1:optimal_value_mark) - optimal_value(optimal_value_mark);
% %if (update_selector ==  1)
%     plot(1:(optimal_value_mark1-1), log(value_diff1(1:optimal_value_mark1-1)), '-*');
%     hold on;
% %else
%     plot(1:(optimal_value_mark2-1), log(value_diff2(1:optimal_value_mark2-1)), '-o');
% %end
% %plot(1:(optimal_value_mark-1), log(value_diff(1:optimal_value_mark-1)), '-*');
% xlabel('Iterations');
% ylabel('ln(F(x) - F(x*))');
% title('Change of ln(F(x) - F(x^*)) w.r.t. the number of iterations');


%% Before figure, let's output the value first to get some sense
%scatter(1:10,iter_plot);
%xlabel('-log_{10}(eps)')
%ylabel('Total number of steps')
%title_str = strcat('m=', num2str(n), ', n=', num2str(p), ', eps=', num2str(input.EPI));
%title(title_str)