clear;clc;close all



n = 100; % 5000 10000 % sample size
K = 100; % 50 100
ng = K; % number of groups
p= 500;

% alpha = 0.1/K * ones(K,1);
% gamma = K / 5; % /20
maxiter = 2000;
gsize=zeros(K,1);
X = randn(n, p);        % A
ind = 1:p;
wtrue = ((-1).^ind) .* exp((1-ind)/100);
wtrue = wtrue(:);
y = X*wtrue + randn(n,1);   % b = Ax+noise
    % Tx=ones(ng,1);  % uniform weight
     
input.A=X;
input.b=y;

%tau=0.1*norm(A'*b,Inf);


L = ng;
lamda = L/5;

% XX = X*X';
% L0 = abs(eigs(XX,1))/gamma/K;
% eps = 2/2/L0;
% mu = min(1/L0, 2*eps);
% Xy = X'*y;              % A'*b
 input.EPI = 1e-4;
input.n = n;
input.p = p;
input.ng = ng;
% report average report
% s_time1=zeros(5,1);
% s_time2=zeros(5,1);
% s_iter1=zeros(5,1);
% s_iter2 = zeros(5,1);
% s_opt1 = zeros(5,1);
% s_opt2 = zeros(5,1);
% s_v1 = zeros(5,1);
% s_v2 = zeros(5,1);
%for i = 1:5
T = randi(2,K,p);  %% Each element is a random number of 1 or 2.
 T = T - 1;  %%  Make it 0 or 1.
 T = sparse(T);
 for j = 1:K
    gsize(j) = length(find(T(j,:) == 1));   % Corrected error: not sum, should be length.
 end
 total_group = sum(gsize);

input.g_size = gsize;
input.T = T;

%[groups] = gen_groupsr(K, T,p);

%fun = @(x)objfunc(x, X, y, lamda, L, lamda1, T);
%x0 = ones(p,1);
%[x, fval] = fminunc(fun, x0);
%fval;

%% Cheng: You need to prepare the accuracy vector
accuracy = zeros(length(1:-0.01:0.01), 1);
i = 0;
for dg = 1:-0.01:(0.01)
    
lamda1 = dg/K;
input.lambda1 = lamda1;
%history5 = rpdmm_overlap_for_random(X,y,K+1,.1,K,groups,gsize,T,alpha,gamma);
%s_opt1 = history5.objval(history5.iter1);
%value_diff1 = abs(history5.objval(1:history5.iter1)-fval);
%s_time1(i) = history5.time1;
%s_iter1(i) = history5.iter1;
% s_opt1 = history5.objval(history5.iter1);
% s_v1=abs((s_opt1(i) - fval)/fval);
%history6 = sadmm_overlap_random(X,y,K+1,.1,K,groups,gsize,T,alpha,gamma);
%history{8} = rpdmm_overlap_for(X,y,1,.1,K,groups,gsize,overlap,alpha,gamma);
[time2,x2,iter_num2,optimal_value_mark2,optimal_value2] = main_alin(input, 2);
%[time2,x2,iter_num2,optimal_value_mark2,optimal_value2] = main_alin(input, 2,s_opt1);
%%value_diff2 = optimal_value2(1:optimal_value_mark2) - fval;
%value_diff2 = abs(optimal_value2(1:iter_num2) - fval);
%s_time2(i) = time2;
%s_iter2(i) = iter_num2;
%  s_opt2 = optimal_value2(optimal_value_mark2);
% s_v2=abs((s_opt2 - fval)/fval);
%plot iterations vs object value
%end
i = i + 1;
accuracy(i) = norm(x2-wtrue);
end
% aver_runtime1 = mean(s_time1);
% aver_runtime2 = mean(s_time2);
% std_r1 = std(s_time1);
% std_r2 = std(s_time2);
% aver_iter1 = mean(s_iter1);
% aver_iter2 = mean(s_iter2);
% stdi1 = std(s_iter1);
% stdi2 = std(s_iter2);
% aver_opt1 = mean(s_opt1);
% aver_opt2 = mean(s_opt2);
% aver_v1 = mean(s_v1);
% aver_v2 = mean(s_v2);
% stdv1 = std(s_v1);
% stdv2 = std(s_v2);
%fprintf('slin: %.2f,%.2f,%.2f,%.2f,%.4f\n', aver_iter2,stdi2,aver_runtime2,std_r2);
%fprintf('pdmm: %.2f,%.2f,%.2f,%.2f,%.4f\n', aver_iter1,stdi1,aver_runtime1,std_r1);

%% Pick the best dg value
[min_norm, index] = min(accuracy);
min_dg_value = 1 - 0.01*(index-1);