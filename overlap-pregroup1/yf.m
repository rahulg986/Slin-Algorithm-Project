function y = yf(x, A, b,L,lamda)  
y = (1/(2*L*lamda))*(norm(A*x-b))^2;
end
