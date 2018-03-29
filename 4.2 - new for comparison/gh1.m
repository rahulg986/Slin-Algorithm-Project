function g_h1 = gh1(x, lambda1) 
g_h1 = x;
g_h1(x>=0) = 1;
g_h1(x<0) = -1;
g_h1 = lambda1*g_h1;
end
