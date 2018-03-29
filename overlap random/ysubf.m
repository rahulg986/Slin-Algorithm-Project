function value = ysubf(x, A, b, sum_gh, xk, d, L, lambda)
value = (1/(2*L*lambda))*(norm(A*x-b))^2 + (sum_gh')*x + 1/2*(x-xk)'*d*(x-xk);
%grad = 1/(input.L*input.lamda)*(input.A'*input.A*x -input.A'*input.b) + input.sum_gh + input.d(1,1)*(x - input.xk);
%grad = grad';
end
