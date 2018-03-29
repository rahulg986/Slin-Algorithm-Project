function [y_hi] = yhi(index,x, lambda1,T)


k = x(T(index,:)==1);
y_hi = lambda1*norm(k);
  
 

end