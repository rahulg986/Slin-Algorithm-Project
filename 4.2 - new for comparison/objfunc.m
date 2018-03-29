function y = objfunc(x, A, b, lambda1, R, lambda2)

y = yf(x, A, b) + yh1(x, lambda1) + yh2(x, R, lambda2);

end