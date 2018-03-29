function [beta,group_ind] = TEST_beta(K,p,m)

% We generate a beta vector of length p, this vector is the joint of K overlapping groups.
% We first (1) generate K-1 groups of pre-assigned values according to the exponential function; Values within the group are of similar order.
%    then  (2) generate location index in the beta vector for each group;
%    then  (3) assign pre-assigned values to these indexes; Deal with overlapping;
%    then  (4) some elements of the beta don’t have values due to overlap, we set this as group K, and assign them all small values. I used exp(-30) here.
%  The final parameters are stored in the vector “beta”.

% K is the number of groups, p is the size of beta. m is the number of elements for each of the first K-1 groups.
% p size of final beta
% K = 8;
% p = 600;
% m is the number of elements per group  m = 100
% so total pre-assign number of elements is m*(K-1)


% step (1)

pre_assign_beta = zeros(K-1,m);
for i = 1 : K-1
    pre_assign_beta(i,:) = ((-1).^((i-1)*m+1:m*i)) .* exp((1-((i-1)*m+1:m*i))/500);
end

% step (2)

group_ind = cell(K,1);
for i = 1:K-1
    group_ind{i} = randsample(p,m)'; % pick 100 index from 1:p range
end

% step (3) 
beta = zeros(1,p);
for i = 1:K-1
    beta(group_ind{i})= pre_assign_beta(i,:);
% this is to assign values to group i locations. If there’s an overlap, new value replaces
% the old one as they are smaller.
end

% step (4)
%find all the left over 0s and their index. Assign them to a small value, and define them as the last group.
last_group_index = find(beta == 0);
beta(last_group_index) = exp(-30);
group_ind{K} = last_group_index;

end


