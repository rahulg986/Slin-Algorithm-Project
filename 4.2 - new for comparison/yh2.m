function y = yh2(x, R, lambda2)

%h2 = lamda_2 *sum_{j} norm(x_{j+1} - x_{j},1)

% y = 0;
% 
% for i = 1:(length(x)-1)
%     y = y + abs(x(i+1) - x(i));
% end
% y = lambda2*y;

y = norm(R*x,1)*lambda2;

end
