clear;
clc;

% Prepare the data
%% 要做的事情：选合适的n和p来测试，n和p决定了矩阵的维度
n=1000;
p=30;
input.n = n;
input.p = p;
x=zeros(p,1);
pos=unidrnd(p,160,1);
for i=1:numel(pos)
    if (mod(pos(i),2)==1) x(pos(i))=1;
    else
        x(pos(i))=-1;
    end
end
input.A=rand(n,p);
A = input.A;
input.b=A*x + normrnd(0,10^(-4),size(A,1),1);
b = input.b
tau=0.1*norm(A'*b,Inf);
input.lambda1 = tau;
input.lambda2 = tau;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 1000;p = 300;
% input.n = n;
% input.p = p;
% input.A = rand(n,p);
% %input.A = input.A(:,1:p);
% input.b = ones(n,1);
% %input.b = input.b(:,1);
% input.lambda1 = 1;
% input.lambda2 = 1;
epi_plot=zeros(10,1);
iter_plot=zeros(10,1);
for i = 1 : 10
    input.EPI = 10^(-i);
    epi_plot(i,1) = input.EPI;
    [x,iter_num] = main_alin(input, 2); 
    iter_plot(i,1)= iter_num;
end
scatter(1:10,iter_plot);
xlabel('-log_{10}(eps)')
ylabel('Total number of steps')
title_str = strcat('m=', num2str(n), ', n=', num2str(p), ', eps=', num2str(input.EPI));
title(title_str)

%input.EPI = 1e-10;  % 根据刚才说的，error的tolerance会根据矩阵的不同而不同，这里设置为可以灵活改变

%x = main_alin(input, 0);
%x = main_alin(input, 2);  % 注意这里只能用1或者2，0是旧的