clear;clc;close all


%overlap=10;  % number of variables overlapped between two consecutive groups     
%p=ng*(g_size-overlap)+overlap;   % totoal number of variables 

n =500; % 5000 10000 % sample size
K =100; % 50 100
ng = K; % number of groups
%gsize = 100;
%overlap = 10;
%d = K*(gsize-overlap) + overlap;
p= 600;
[wtrue,group_ind]=TEST_beta(100,600,150);
alpha = 0.01/K * ones(K,1);
gamma = K / 5; % /20
maxiter = 2000;
gsize=zeros(K,1);
s_time1=zeros(10,1);
s_time2=zeros(10,1);
s_iter1=zeros(10,1);
s_iter2 = zeros(10,1);

for l = 1:10
X = randn(n, p);        % A
ind = 1:p;
%wtrue = ((-1).^ind) .* exp((1-ind)/100);
wtrue = wtrue(:);
y = X*wtrue + randn(n,1);   % b = Ax+noise
    % Tx=ones(ng,1);  % uniform weight
     
input.A=X;
input.b=y;
input.group_ind = group_ind;
%tau=0.1*norm(A'*b,Inf);
lamda1 = 0.01/K;
input.lambda1 = lamda1;

L = ng;
lamda = L/5;

XX = X*X';
L0 = abs(eigs(XX,1))/gamma/K;
eps = 2/2/L0;
mu = min(1/L0, 2*eps);
Xy = X'*y;              % A'*b
 input.EPI = 1e-3;
input.n = n;
input.p = p;
input.ng = ng;
% report average report

% s_opt1 = zeros(10,1);
% s_opt2 = zeros(10,1);
% s_v1 = zeros(10,1);
% s_v2 = zeros(10,1);

%T = randi(2,K,p);  %% Each element is a random number of 1 or 2.
% T = T - 1;  %%  Make it 0 or 1.
% T = sparse(T);
%% New thing: get the T from the group_ind
T = zeros(K, p);
for j = 1:K
    for i = 1:length(group_ind{j})
        T(j,group_ind{j}(i)) = 1;
    end
end

 for j = 1:K
    gsize(j) = length(find(T(j,:) == 1));   % Corrected error: not sum, should be length.
 end
 total_group = sum(gsize);

input.g_size = gsize;
input.T = T;

%[groups] = gen_groups(K, gsize, overlap);
[groups] = gen_groupsr(K, T,p);

%fun = @(x)objfunc(x, X, y, lamda, L, lamda1, T);
%x0 = ones(p,1);
%[x, fval] = fminunc(fun, x0);
%fval;
% x = randn(p, 1);


%history{5} = rpdmm_overlap_for(X,y,K+1,.1,K,groups,gsize,overlap,alpha,gamma);

history5 = rpdmm_overlap_for_random(X,y,K+1,.1,K,groups,gsize,T,alpha,gamma);
s_opt1 = history5.objval(history5.iter1);
%value_diff1 = abs(history5.objval(1:history5.iter1)-fval);
s_time1(l) = history5.time1;
s_iter1(l) = history5.iter1;
% s_opt1 = history5.objval(history5.iter1);
% s_v1=abs((s_opt1(i) - fval)/fval);
%history6 = sadmm_overlap_random(X,y,K+1,.1,K,groups,gsize,T,alpha,gamma);
%history{8} = rpdmm_overlap_for(X,y,1,.1,K,groups,gsize,overlap,alpha,gamma);
%[time2,x2,iter_num2,optimal_value_mark2,optimal_value2] = main_alin(input, 2);
[time2,x2,iter_num2,optimal_value_mark2,optimal_value2] = main_alin(input, 2,s_opt1);
%%value_diff2 = optimal_value2(1:optimal_value_mark2) - fval;
%value_diff2 = abs(optimal_value2(1:iter_num2) - fval);
s_time2(l) = time2;
s_iter2(l) = iter_num2;
%  s_opt2 = optimal_value2(optimal_value_mark2);
% s_v2=abs((s_opt2 - fval)/fval);
%plot iterations vs object value
end
aver_runtime1 = mean(s_time1);
aver_runtime2 = mean(s_time2);
std_r1 = std(s_time1);
std_r2 = std(s_time2);
aver_iter1 = mean(s_iter1);
aver_iter2 = mean(s_iter2);
stdi1 = std(s_iter1);
stdi2 = std(s_iter2);
% aver_opt1 = mean(s_opt1);
% aver_opt2 = mean(s_opt2);
% aver_v1 = mean(s_v1);
% aver_v2 = mean(s_v2);
% stdv1 = std(s_v1);
% stdv2 = std(s_v2);
fprintf('slin: %.4f,%.4f,%.4f,%.4f,%.4f\n', aver_iter2,stdi2,aver_runtime2,std_r2);
fprintf('pdmm: %.4f,%.4f,%.4f,%.4f,%.4f\n', aver_iter1,stdi1,aver_runtime1,std_r1);



% %% Call fminunc
% fun = @(x)objfunc(x, X, y, lamda, L, lamda1, T);
% x0 = ones(p,1);
% [x, fval] = fminunc(fun, x0);
% fval