clear;clc;close all

n = 5000; % 5000 10000
d = 1000;

[X, y, wtrue] = gentoy_graph(n, d);
opts = struct('cortype', 1, 'corthreshold', 0.7);
[C, CNorm, E] = gennetwork(X, opts);
K = size(E,1);
alpha = 1/K * ones(K,1);
gamma = 100;
maxiter = 1000;

XX = X'*X;
L0 = abs(eigs(XX,1))/gamma/K;
eps = 0.0001;
mu = min(1/L0, 2*eps);
Xy = X'*y;

x = randn(d, 1);
aux_pa = @(x) aux_graph(x, XX, Xy, y, 1/gamma/K, alpha, mu, E);
% aux_s = @(x) aux_graph(x, XX, Xy, y, 1/gamma/K, alpha, 2*eps, E);

% [w_pa_pg, f_pa_pg] = pa_pg(x, mu, aux_pa, maxiter);
% [w_s_pg, f_s_pg] = s_pg(x, 2*eps, L0, aux_s, maxiter);

[w_pa_apg, f_pa_apg] = pa_apg(x, mu, aux_pa, maxiter);
% [w_s_apg, f_s_apg] = s_apg(x, 2*eps, L0, aux_s, maxiter);


figure;
% semilogy(1:maxiter, f_pa_pg, ':m', 'LineWidth',2); hold on;
% semilogy(1:maxiter, f_s_pg, '-.c', 'linewidth', 2); hold on;
semilogy(1:maxiter, f_pa_apg, 'r', 'linewidth', 2); hold on;
% semilogy(1:maxiter, f_s_apg, '--b', 'linewidth', 2); hold on;
% legend('PA-PG', 'S-PG', 'PA-APG', 'S-APG', 'fontsize', 20)
% %     set(findall(fhandle,'type','text'),'fontSize',14,'fontWeight','bold');
% title('2\epsilon = 2/L_0', 'fontsize', 20);
% set(0,'DefaultAxesFontSize',20)