function y = min_objfunc(x, z_f, z_h1, z_h2, A, b, lambda1, lambda2, g_f, g_h1, g_h2)
value1 = yh1(x, lambda1) - (yh1(z_h1, lambda1) + (g_h1')*(x - z_h1));
value2 = yh2(x, lambda2) - (yh2(z_h2, lambda2) + (g_h2')*(x - z_h2));
if (value1> value2)
     y = yf(z_f, A, b) + (g_f')*(x-z_f) + yh1(z_h1, lambda1) + (g_h1')*(x - z_h1) + yh2(x, lambda2);
else
   y = yf(z_f, A, b) + (g_f')*(x-z_f) + yh1(x, lambda1) + yh2(z_h2, lambda2) + (g_h2')*(x - z_h2);
end
end