function [value, grad] = ysubf(x, input)
x = x';
value = (1/(2*input.L*input.lamda))*(norm(input.A*x-input.b))^2 + input.sum_gh'*x + 1/2*input.d(1,1)*norm(x - input.xk)^2;
grad = 1/(input.L*input.lamda)*(input.A'*input.A*x -input.A'*input.b) + input.sum_gh + input.d(1,1)*(x - input.xk);
grad = grad';
end
