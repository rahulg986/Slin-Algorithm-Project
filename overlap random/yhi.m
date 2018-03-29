%function [y_hi] = yhi(index,x, lambda1,T,group_ind)
function [y_hi] = yhi(index,x, lambda1,group_ind)

%k = x(T(index,:)==1);
%y_hi = lambda1*norm(k);
 % k=group_ind{index};
 y_hi = lambda1*norm(x(group_ind{index}));

end