function g_h2 = gh2(x, R, lambda2)

% n = length(x);
% g_h2 = zeros(n,1);
% g_h2(1) = sign(x(1) - x(2));
% for i = 2:(n-1)
%     g_h2(i) = sign(x(i) - x(i-1)) + sign(x(i) - x(i+1));
% end
% g_h2(n) = sign(x(n)-x(n-1));
% 
% g_h2 = lambda2*g_h2;

y = R*x;
p = length(x);
g_h2 = zeros(p,1);
k = length(y);
for i = 1:p
    for j = 1:k
        if (y(j) >= 0)
            g_h2(i) = g_h2(i) + R(j,i);
        else
            g_h2(i) = g_h2(i) - R(j,i);
        end
    end
end

g_h2 = g_h2*lambda2;
end