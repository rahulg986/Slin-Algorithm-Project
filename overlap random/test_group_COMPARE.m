clear;clc;close all


%overlap=10;  % number of variables overlapped between two consecutive groups     
%p=ng*(g_size-overlap)+overlap;   % totoal number of variables 

n =500; % 5000 10000 % sample size
K = 200; % 50 100
ng = K; % number of groups
%gsize = 100;
%overlap = 10;
%d = K*(gsize-overlap) + overlap;
p= 600;
alpha = 0.1/K * ones(K,1);
gamma = K / 5; % /20
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
lamda1 = 0.1/K;
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
s_time1=zeros(10,1);
s_time2=zeros(10,1);
s_iter1=zeros(10,1);
s_iter2 = zeros(10,1);
% s_opt1 = zeros(10,1);
% s_opt2 = zeros(10,1);
% s_v1 = zeros(10,1);
% s_v2 = zeros(10,1);
for i = 1:10
T = randi(2,K,p);  %% Each element is a random number of 1 or 2.
 T = T - 1;  %%  Make it 0 or 1.
 T = sparse(T);
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
%  aux_pa = @(x) aux_group(x, X, Xy, y, 1/gamma/K, alpha, mu, groups);
%  aux_s = @(x) aux_group(x, X, Xy, y, 1/gamma/K, alpha, 2*eps, groups);
% % 
% history{1} = pa_pg(x, mu, aux_pa, maxiter);
 
% history{2} = s_pg(x, 2*eps, L0, aux_s, maxiter);
% 
% history{3} = pa_apg(x, mu, aux_pa, maxiter);
% 
% history{4} = s_apg(x, 2*eps, L0, aux_s, maxiter);

%history{5} = rpdmm_overlap_for(X,y,K+1,.1,K,groups,gsize,overlap,alpha,gamma);

history5 = rpdmm_overlap_for_random(X,y,K+1,.1,K,groups,gsize,T,alpha,gamma);
s_opt1 = history5.objval(history5.iter1);
%value_diff1 = abs(history5.objval(1:history5.iter1)-fval);
s_time1(i) = history5.time1;
s_iter1(i) = history5.iter1;
% s_opt1 = history5.objval(history5.iter1);
% s_v1=abs((s_opt1(i) - fval)/fval);
%history6 = sadmm_overlap_random(X,y,K+1,.1,K,groups,gsize,T,alpha,gamma);
%history{8} = rpdmm_overlap_for(X,y,1,.1,K,groups,gsize,overlap,alpha,gamma);
%[time2,x2,iter_num2,optimal_value_mark2,optimal_value2] = main_alin(input, 2);
[time2,x2,iter_num2,optimal_value_mark2,optimal_value2] = main_alin(input, 2,s_opt1);
%%value_diff2 = optimal_value2(1:optimal_value_mark2) - fval;
%value_diff2 = abs(optimal_value2(1:iter_num2) - fval);
s_time2(i) = time2;
s_iter2(i) = iter_num2;
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
fprintf('slin: %.2f,%.2f,%.2f,%.2f,%.4f\n', aver_iter2,stdi2,aver_runtime2,std_r2);
fprintf('pdmm: %.2f,%.2f,%.2f,%.2f,%.4f\n', aver_iter1,stdi1,aver_runtime1,std_r1);
% %plot iteration vs ln|f -f*|
%  plot(1:(iter_num2-1), log(value_diff2(1:iter_num2-1)), '--');
%  ylabel('ln(F(x) - F(x*))');
%  xlabel('Iterations');
%  hold on
%  plot(1:length(history5.objval), log(value_diff1(1:length(history5.objval))),'-');
%  ylabel('ln(F(x) - F(x*))');
%  xlabel('Iterations');

% %plot iteration vs object value
% semilogy(1:optimal_value_mark2,optimal_value2(1:optimal_value_mark2),'-');
% ylabel('Objective values');
% xlabel('Iterations');
% hold on
% semilogy(1:length(history5.objval), history5.objval,'--');
% ylabel('Objective values');
% xlabel('Iterations');
% %plot time vs object value
% hold on
% semilogy(time2,optimal_value2(1:iter_num2),'-');
% ylabel('Objective values');
% xlabel('Iterations');
% hold on
% semilogy(history5.time, history5.objval,'--');
% ylabel('Objective values');
% xlabel('Iterations');
%history{9} = rpdmm_overlap_for(X,y,11,.1,K,groups,gsize,overlap,alpha,gamma);

%history{10} = rpdmm_overlap_for(X,y,21,.1,K,groups,gsize,overlap,alpha,gamma);

%history{11} = rpdmm_overlap_for(X,y,31,.1,K,groups,gsize,overlap,alpha,gamma);

%history{12} = rpdmm_overlap_for(X,y,41,.1,K,groups,gsize,overlap,alpha,gamma);

%history{13} = rpdmm_overlap_for(X,y,51,.1,K,groups,gsize,overlap,alpha,gamma);

%history{14} = rpdmm_overlap_for(X,y,61,.1,K,groups,gsize,overlap,alpha,gamma);

%history{15} = rpdmm_overlap_for(X,y,71,.1,K,groups,gsize,overlap,alpha,gamma);

%history{16} = rpdmm_overlap_for(X,y,81,.1,K,groups,gsize,overlap,alpha,gamma);%

%history{17} = rpdmm_overlap_for(X,y,91,.1,K,groups,gsize,overlap,alpha,gamma);

% history{6} = admm_overlap_for(X,y,K+1,.1,K,groups,gsize,overlap,alpha,gamma);
% 
% history{7} = sadmm_overlap(X,y,K+1,.1,K,groups,alpha,gamma);

%history{18} = rbsumm_overlap_for(X,y,.1,K,groups,gsize,overlap,alpha,gamma);

% figure; hold on
% for i = 1:7
%     plot(history{i}.time, history{i}.objval);
%     plot(1:length(history{i}.time), log(history{i}.objval));
% end
% 
% figure; hold on
% plot(history{5}.time, log10(history{5}.r_norm + history{5}.s_norm));
% plot(1:length(history{5}.time), log10(history{5}.r_norm + history{5}.s_norm));
% for i = 8:2:17
%     plot(history{i}.time, log10(history{i}.r_norm + history{i}.s_norm));
%     plot(1:length(history{i}.time), log10(history{i}.r_norm + history{i}.s_norm));
% end
% 
% figure; hold on
% for i = 1:7
%     plot(1:length(history{i}.time), history{i}.objval);
% end
% 
% figure; hold on
% for i = 1:4
%     plot(history{i}.time, log10(history{i}.res));
% end
% 
% for i = 5:7
%     plot(history{i}.time, log10(history{i}.r_norm + history{i}.s_norm));
% end
% 
% figure; hold on
% for i = 1:4
%     plot(1:length(history{i}.time), log10(history{i}.res));
% end
% 
% for i = 5:7
%     plot(1:length(history{i}.time), log10(history{i}.r_norm + history{i}.s_norm));
% end

% fs = [f_pa_pg(:), f_pa_apg(:), f_s_pg(:), f_s_apg(:)];
% figure;
% semilogy(1:maxiter, fs(1:end,:));

% figure;
% semilogy(1:maxiter, f_pa_pg, ':m', 'LineWidth',2); hold on;
% semilogy(1:maxiter, f_s_pg, '-.c', 'linewidth', 2); hold on;
% semilogy(1:maxiter, f_pa_apg, 'r', 'linewidth', 2); hold on;
% semilogy(1:maxiter, f_s_apg, '--b', 'linewidth', 2); hold on;
% legend('PA-PG', 'S-PG', 'PA-APG', 'S-APG', 'fontsize', 20)
% %     set(findall(fhandle,'type','text'),'fontSize',14,'fontWeight','bold');
% title('2\epsilon = 2/L_0', 'fontsize', 20);
% set(0,'DefaultAxesFontSize',20)

% %% Call fminunc
% fun = @(x)objfunc(x, X, y, lamda, L, lamda1, T);
% x0 = ones(p,1);
% [x, fval] = fminunc(fun, x0);
% fval