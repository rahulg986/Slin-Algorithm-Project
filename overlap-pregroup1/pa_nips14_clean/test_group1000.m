clear;clc;close all

n = 5000; % 5000 10000
K = 10000; % 50 100
gsize = 50;
overlap = 10;
d = K*(gsize-overlap) + overlap;

alpha = 1/K * ones(K,1);
gamma = K / 5; % /20
maxiter = 2000;

X = randn(n, d);        % A
ind = 1:d;
wtrue = ((-1).^ind) .* exp((1-ind)/100);
wtrue = wtrue(:);
y = X*wtrue + randn(n,1);   % b = Ax+noise

XX = X*X';
L0 = abs(eigs(XX,1))/gamma/K;
eps = 2/2/L0;
mu = min(1/L0, 2*eps);
Xy = X'*y;              % A'*b

[groups] = gen_groups(K, gsize, overlap);

x = randn(d, 1);
aux_pa = @(x) aux_group(x, X, Xy, y, 1/gamma/K, alpha, mu, groups);
aux_s = @(x) aux_group(x, X, Xy, y, 1/gamma/K, alpha, 2*eps, groups);

history = pdmm_overlap_for(X,y,K+1,.1,K,groups,gsize,overlap,alpha,gamma);

history = rpdmm_overlap_for(X,y,K+1,.1,K,groups,gsize,overlap,alpha,gamma);

history = rpdmm_overlap_for(X,y,1,.1,K,groups,gsize,overlap,alpha,gamma);

history = rpdmm_overlap_for(X,y,101,.1,K,groups,gsize,overlap,alpha,gamma);

history = rpdmm_overlap_for(X,y,201,.1,K,groups,gsize,overlap,alpha,gamma);

history = rpdmm_overlap_for(X,y,301,.1,K,groups,gsize,overlap,alpha,gamma);

history = rpdmm_overlap_for(X,y,401,.1,K,groups,gsize,overlap,alpha,gamma);

history = rpdmm_overlap_for(X,y,501,.1,K,groups,gsize,overlap,alpha,gamma);

history = rpdmm_overlap_for(X,y,601,.1,K,groups,gsize,overlap,alpha,gamma);

history = rpdmm_overlap_for(X,y,701,.1,K,groups,gsize,overlap,alpha,gamma);

history = rpdmm_overlap_for(X,y,801,.1,K,groups,gsize,overlap,alpha,gamma);

history = rpdmm_overlap_for(X,y,901,.1,K,groups,gsize,overlap,alpha,gamma);

history = rpdmm_overlap_for(X,y,1001,.1,K,groups,gsize,overlap,alpha,gamma);