function [groups, A] = gen_groups(K, gsize, overlap)

groups = cell(K,1);
groups{1} = 1:gsize;
idx = gsize - overlap + 1;
A = zeros(K, K*(gsize-overlap)+overlap);
for ii = 2:K
    groups{ii} = idx:(idx+gsize-1);
    idx = idx + gsize - overlap;
end