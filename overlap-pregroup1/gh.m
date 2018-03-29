function g_h = gh(x, lambda1,T) 
p = length(T);
g_h = zeros(p,1);

m=sumsqr(x(T==1));

for j = 1:p   
    if (T(j)==0)
        g_h(j)=0;
    else
        g_h(j) = x(j)/sqrt(m);
    end
end

g_h = lambda1*g_h;
end
