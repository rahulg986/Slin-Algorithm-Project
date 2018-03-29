function [groups, A] = gen_groupsr(K, T,d)

groups = cell(K,1);
for ii=1:K
    groups{ii} = find(T(ii,:));
end

%groups{1} = 1:gsize;
%idx = gsize - overlap + 1;
%A = zeros(K, K*(gsize-overlap)+overlap);
A = zeros(K, d);
% for ii = 2:K
%     groups{ii} = idx:(idx+gsize-1);
%     idx = idx + gsize - overlap;
% end