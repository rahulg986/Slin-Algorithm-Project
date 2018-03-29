clear;
clc;

%% Prepare the data
%n = 3000;
n_list=[100]; %list of n values to try
p=10;

%p_list = [500,1000,3000,4000];  % List of p values to try
% time array for each method
time2list= zeros(length(n_list),1);
time1list= zeros(length(n_list),1);
time3list= zeros(length(n_list),1);
time_1list= zeros(length(n_list),1);
time_2list= zeros(length(n_list),1);
time_3list= zeros(length(n_list),1);
time_4list= zeros(length(n_list),1);

% time2list= zeros(length(p_list),1);
% time1list= zeros(length(p_list),1);
% time3list= zeros(length(p_list),1);      
% time_1list= zeros(length(p_list),1);
% time_2list= zeros(length(p_list),1);
% time_3list= zeros(length(p_list),1);
% time_4list= zeros(length(p_list),1);
for j = 1:length(n_list) %go through all n values
%for j = 1:length(p_list)  % go through all p values
%p= p_list(j);
n = n_list(j);
input.EPI = 1e-3;
input.n = n;
input.p = p;
x=zeros(p,1);
pos=unidrnd(p,160,1);
for i=1:numel(pos)
    if (mod(pos(i),2)==1) 
        x(pos(i))=1;
    else
        x(pos(i))=-1;
    end
end
R = zeros(p-1,p);
for i = 1:(p-1)
    R(i,i) = -1;
    R(i,i+1) = 1;
end

input.A=rand(n,p);
%for sparsity
for k = 1:p
    s= randperm(n);
s= s(1:0.8*n);
   input.A(s,k)=0;
end
%for normalize A
% for k= 1:p
%    input.A(:,k)= C(:,k)/norm(C(:,k));
% end
A = input.A;
input.b=A*x + normrnd(0,10^(-4),size(A,1),1);
b = input.b;
tau=0.1*norm(A'*b,Inf);
input.lambda1 =  0.5;%0.01*tau
input.lambda2 = 0.5;
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
%%%%%%%%
% fun = @(x)objfunc(x, A, b, input.lambda1, R, input.lambda2);
% x0 = ones(p,1);
% [x, fval] = fminunc(fun, x0);
% fval;
%% Call the four methods for comparison

[time2list(j), v2, v2_m, x2, iter_num2, optimal_value_mark2, optimal_value2] = main_alin(input, 2);%slin

[time1list(j), v1, v1_m, x1, iter_num1, optimal_value_mark1, optimal_value1] = main_alin(input, 1);%alin
%[time0, v0,v0_m,x0,iter_num0,optimal_value_mark0,optimal_value0] = main_alin(input, 0); %peaceman
%time2list(j) = time2;
[time_1list(j), v_1, v_1_m, x_1, iter_num_1, optimal_value_mark_1, optimal_value_1] = main_alin(input, -1); %douglas

[time_2list(j), v_2, v_2_m, x_2, iter_num_2, optimal_value_mark_2, optimal_value_2] = main_alin(input, -2); %douglas over relax
[time_3list(j), v_3, v_3_m, x_3, iter_num_3, optimal_value_mark_3, optimal_value_3] = main_alin(input, -3); %douglas under relax
[time_4list(j), v_4, v_4_m, x_4, iter_num_4, optimal_value_mark_4, optimal_value_4] = main_alin(input, -4); %douglas with \rho

[time3list(j), x3, iter_num3, optimal_value_mark3, optimal_value3] = bangcongvu(input);
%[time3list(j), x3, iter_num3, optimal_value_mark3, optimal_value3] = bangcongvu(input, time2list(j) * 2); %call bangcongvu
end

%  fprintf('selective: %.2f,%.2f,%.2f,%.2f,%.4f,%.4f\n', aver_iter2,stdi2,aver_runtime2,std_r2,aver_v2,stdv2);
%   fprintf('bundle: %.2f,%.2f,%.2f,%.2f,%.4f,%.4f\n', aver_iter1,stdi1,aver_runtime1,std_r1,aver_v1,stdv1);
%    fprintf('peaceman: %.2f,%.2f,%.2f,%.2f,%.4f,%.4f\n', aver_iter0,stdi0,aver_runtime0,std_r0,aver_v0,stdv0);
%  fprintf('douglas: %.2f,%.2f,%.2f,%.2f,%.4f,%.4f\n', aver_iter_1,stdi_1,aver_runtime_1,std_r_1,aver_v_1,stdv_1);
%  fprintf('douglasrelax: %.2f,%.2f,%.2f,%.2f,%.4f,%.4f\n', aver_iter_2,stdi_2,aver_runtime_2,std_r_2,aver_v_2,stdv_2);
%  fprintf('bangcongvu: %.2f,%.2f,%.2f,%.2f,%.4f,%.4f\n', aver_iter3,stdi3,aver_runtime3,std_r3,aver_v3,stdv3);


%% Plot computational time over dimension of variables
% figure;
% plot(p_list, time2list);
% hold on
% plot(p_list, time1list);
% hold on
% plot(p_list, time_1list);
% hold on
% plot(p_list, time_2list);
% hold on
% plot(p_list, time_3list);
% hold on
% plot(p_list, time_4list);
% hold on
% plot(p_list, time3list);
% ylabel('Running time');
% xlabel('Variable dimension');
% title('Change of running time w.r.t. the dimension of variables');
% hold off

%% Plot computational time over size of the problem
% figure;
% plot(n_list, time2list);
% hold on
% plot(n_list, time1list);
% hold on
% plot(n_list, time_1list);
% hold on
% plot(n_list, time_2list);
% hold on
% plot(n_list, time_3list);
% hold on
% plot(n_list, time_4list);
% hold on
% plot(n_list, time3list);
% ylabel('Running time');
% xlabel('Sample size');
% title('Change of running time w.r.t. the sample size');
% hold off
% 
% % %% plot for object value vs iteration and time
% %  figure
% % % %Plot the optimal value changes with respect to the number of iterations
% %  loglog(1:optimal_value_mark2,optimal_value2(1:optimal_value_mark2), '-k');
% %  hold on
% %  loglog(1:optimal_value_mark1,optimal_value1(1:optimal_value_mark1),'-b');
% %  hold on
% %  loglog(1:optimal_value_mark_1,optimal_value_1(1:optimal_value_mark_1),'-g');
% %  hold on
% %  loglog(1:optimal_value_mark_2,optimal_value_2(1:optimal_value_mark_2),'-y');
% %  hold on
% %  loglog(1:optimal_value_mark3,optimal_value3(1:optimal_value_mark3),'-r');
% %  ylabel('objective value');
% %  xlabel('Iterations');
% %  title('Change of objective value w.r.t. the number of iterations');
%   
% % %figure ln(f(x)-f(x*)) with slin as benchmark
value_diff2= optimal_value2(1:optimal_value_mark2) - optimal_value2(optimal_value_mark2);
value_diff1= optimal_value1(1:optimal_value_mark1) - optimal_value2(optimal_value_mark2);
value_diff3= optimal_value3(1:optimal_value_mark3) - optimal_value2(optimal_value_mark2);
value_diff_1= optimal_value_1(1:optimal_value_mark_1) - optimal_value2(optimal_value_mark2);
value_diff_2= optimal_value_2(1:optimal_value_mark_2) - optimal_value2(optimal_value_mark2);
value_diff_3= optimal_value_3(1:optimal_value_mark_3) - optimal_value2(optimal_value_mark2);
value_diff_4= optimal_value_4(1:optimal_value_mark_4) - optimal_value2(optimal_value_mark2);
figure
semilogx(1:(optimal_value_mark2-1), log(value_diff2(1:optimal_value_mark2-1)));
hold on
semilogx(1:(optimal_value_mark1-1), log(value_diff1(1:optimal_value_mark1-1)));
hold on
semilogx(1:(optimal_value_mark_1-1), log(value_diff_1(1:optimal_value_mark_1-1)));
hold on
semilogx(1:(optimal_value_mark_2-1), log(value_diff_2(1:optimal_value_mark_2-1)));
hold on
semilogx(1:(optimal_value_mark_3-1), log(value_diff_3(1:optimal_value_mark_3-1)));
hold on
semilogx(1:(optimal_value_mark_4-1), log(value_diff_4(1:optimal_value_mark_4-1)));
hold on
semilogx(1:(optimal_value_mark3-1), log(value_diff3(1:optimal_value_mark3-1)));
hold on
  xlabel('Iterations');
  ylabel('ln(F(x) - F(x^*))');
  title('Change of ln(F(x) - F(x^*)) w.r.t. the number of iterations');
 hold off
%%%%%%%%figure ln(f(x)-f(x*))
%  figure
% % %value_diff = optimal_value(1:optimal_value_mark) - optimal_value(optimal_value_mark);
%  if (update_selector ==  1)
%      plot(1:(optimal_value_mark1-1), log(value_diff1(1:optimal_value_mark1-1)), '-');
%      hold on;
%  else
%      plot(1:(optimal_value_mark2-1), log(value_diff2(1:optimal_value_mark2-1)), '--');
%  end
% % 
% % %plot(1:(optimal_value_mark-1), log(value_diff(1:optimal_value_mark-1)), '-');
%  xlabel('Iterations');
%  ylabel('ln(F(x) - F(x*))');
%  title('Change of ln(F(x) - F(x^*)) w.r.t. the number of iterations');
% % 

%% Before figure, let's output the value first to get some sense
%scatter(1:10,iter_plot);
%xlabel('-log_{10}(eps)')
%ylabel('Total number of steps')
%title_str = strcat('m=', num2str(n), ', n=', num2str(p), ', eps=', num2str(input.EPI));
%title(title_str)