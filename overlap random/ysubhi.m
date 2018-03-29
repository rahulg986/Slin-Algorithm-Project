function value = ysubhi(x, next_index, a_f, a_h, T, xk, d, lambda1)
l = a_f + sum(a_h,2) - a_h(:, next_index);
T_1 = T(next_index, :) == 1;
value = lambda1*sqrt(sumsqr(x(T_1))) + (l')*x + 1/2*(x-xk)'*d*(x-xk);
end
