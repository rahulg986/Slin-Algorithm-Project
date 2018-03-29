clear;clc;close all


%overlap=10;  % number of variables overlapped between two consecutive groups     
%p=ng*(g_size-overlap)+overlap;   % totoal number of variables 

n = 1000; % 5000 10000 % sample size
K = 100; % 50 100
ng = K; % number of groups
%gsize = 100;
%overlap = 10;
%d = K*(gsize-overlap) + overlap;
d=100;
p = d;
alpha = 1/K * ones(K,1);
gamma = K / 5; % /20
maxiter = 2000;
gsize=zeros(K,1);
T = randi(2,K,d);  %% Each element is a random number of 1 or 2.
 T = T - 1;  %%  Make it 0 or 1.
 T = sparse(T);
 for i = 1:K
    gsize(i) = sum(find(T(i,:)));
 end
 total_group = sum(gsize);
 input.EPI = 1e-3;
input.n = n;
input.p = p;
input.ng = ng;
input.g_size = gsize;

X = randn(n, d);        % A
ind = 1:d;
wtrue = ((-1).^ind) .* exp((1-ind)/100);
wtrue = wtrue(:);
y = X*wtrue + randn(n,1);   % b = Ax+noise
    % Tx=ones(ng,1);  % uniform weight
     
input.A=X;
input.b=y;
input.T = T;
%tau=0.1*norm(A'*b,Inf);
lamda1 = 1/K;
input.lambda1 = lamda1;



XX = X*X';
L0 = abs(eigs(XX,1))/gamma/K;
eps = 2/2/L0;
mu = min(1/L0, 2*eps);
Xy = X'*y;              % A'*b

%[groups] = gen_groups(K, gsize, overlap);
[groups] = gen_groupsr(K, T,d);
x = randn(d, 1);
 aux_pa = @(x) aux_group(x, X, Xy, y, 1/gamma/K, alpha, mu, groups);
 aux_s = @(x) aux_group(x, X, Xy, y, 1/gamma/K, alpha, 2*eps, groups);
% 
% history{1} = pa_pg(x, mu, aux_pa, maxiter);
 
% history{2} = s_pg(x, 2*eps, L0, aux_s, maxiter);
% 
% history{3} = pa_apg(x, mu, aux_pa, maxiter);
% 
% history{4} = s_apg(x, 2*eps, L0, aux_s, maxiter);

%history{5} = rpdmm_overlap_for(X,y,K+1,.1,K,groups,gsize,overlap,alpha,gamma);
history{5} = rpdmm_overlap_for_random(X,y,K+1,.1,K,groups,gsize,T,alpha,gamma);
%history{8} = rpdmm_overlap_for(X,y,1,.1,K,groups,gsize,overlap,alpha,gamma);
[x2,iter_num2,optimal_value_mark2,optimal_value2] = main_alin(input, 2);
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