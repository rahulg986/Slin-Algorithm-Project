function g_f = gf(x, A, b, L, lamda)
g_f = (1/(lamda*L))*(A')*A*x - (A')*b;
end

