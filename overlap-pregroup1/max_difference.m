function [next_index, total_diff] = max_difference(x, a_f,b_f, a_h,b_h,T, A, b, lamda, L,lambda1,next_index)
next_index = 0;
f_diff = yf(x, A, b,L,lamda) - (a_f'*x +b_f);
total_diff = f_diff;
ng = size(T,1); 
temp_diff = f_diff;
for i = 1:ng
    hi_diff= yhi(i,x, lambda1, T) - (a_h(:,i)'*x + b_h(i));
    total_diff = total_diff + hi_diff;
    if (temp_diff < hi_diff)
        temp_diff = hi_diff;
        next_index = i;
    end

end
    
    
