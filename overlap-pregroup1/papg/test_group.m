clear;clc;close all

n = 5000; % 5000 10000
K = 50; % 50 100
gsize = 100;
overlap = 10;
d = K*(gsize-overlap) + overlap;

alpha = 1/K * ones(K,1);
gamma = K / 5; % /20
maxiter = 150;

X = randn(n, d);
ind = 1:d;
wtrue = ((-1).^ind) .* exp((1-ind)/100);
wtrue = wtrue(:);
y = X*wtrue + randn(n,1);

XX = X'*X;
L0 = abs(eigs(XX,1))/gamma/K;
eps = 2/2/L0;
mu = min(1/L0, 2*eps);
Xy = X'*y;

[groups] = gen_groups(K, gsize, overlap);

x = randn(d, 1);
aux_pa = @(x) aux_group(x, XX, Xy, y, 1/gamma/K, alpha, mu, groups);
aux_s = @(x) aux_group(x, XX, Xy, y, 1/gamma/K, alpha, 2*eps, groups);

[w_pa_pg, f_pa_pg] = pa_pg(x, mu, aux_pa, maxiter);
[w_s_pg, f_s_pg] = s_pg(x, 2*eps, L0, aux_s, maxiter);

[w_pa_apg, f_pa_apg] = pa_apg(x, mu, aux_pa, maxiter);
[w_s_apg, f_s_apg] = s_apg(x, 2*eps, L0, aux_s, maxiter);

% fs = [f_pa_pg(:), f_pa_apg(:), f_s_pg(:), f_s_apg(:)];
% figure;
% semilogy(1:maxiter, fs(1:end,:));

figure;
semilogy(1:maxiter, f_pa_pg, ':m', 'LineWidth',2); hold on;
semilogy(1:maxiter, f_s_pg, '-.c', 'linewidth', 2); hold on;
semilogy(1:maxiter, f_pa_apg, 'r', 'linewidth', 2); hold on;
semilogy(1:maxiter, f_s_apg, '--b', 'linewidth', 2); hold on;
legend('PA-PG', 'S-PG', 'PA-APG', 'S-APG', 'fontsize', 20)
%     set(findall(fhandle,'type','text'),'fontSize',14,'fontWeight','bold');
title('2\epsilon = 2/L_0', 'fontsize', 20);
set(0,'DefaultAxesFontSize',20)