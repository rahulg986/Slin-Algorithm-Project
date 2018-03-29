function [y_h] = yh(x, lambda1,T)
ng = size(T,1);
y_h = 0;

for i = 1:ng
    k = x(T(i,:)==1);
    temp_y = lambda1*norm(k);
    y_h = y_h+temp_y;

end

end