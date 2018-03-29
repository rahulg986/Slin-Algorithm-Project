clear;
clc;

%% Prepare the data
n=300;
p=1000;
input.EPI = 1e-7;
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
s_time_1=zeros(5,1);
s_time0=zeros(5,1);
s_time1=zeros(5,1);
s_time2=zeros(5,1);
s_iter_1=zeros(5,1);
s_iter0=zeros(5,1);
s_iter1 = zeros(5,1);
s_iter2 = zeros(5,1);
s_opt_1=zeros(5,1);
s_opt0 = zeros(5,1);
s_opt1 = zeros(5,1);
s_opt2 = zeros(5,1);
s_v_1=zeros(5,1);
s_v0 = zeros(5,1);
s_v1 = zeros(5,1);
s_v2 = zeros(5,1);
for i = 1:5
input.A=rand(n,p);
A = input.A;
input.b=A*x + normrnd(0,10^(-4),size(A,1),1);
b = input.b;
tau=0.1*norm(A'*b,Inf);
input.lambda1 = 0.01*tau;
input.lambda2 = 0.01*tau;
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
fun = @(x)objfunc(x, A, b, input.lambda1, R, input.lambda2);
x0 = ones(p,1);
[x, fval] = fminunc(fun, x0);
fval;
%% Call the four methods for comparison

  [time_1,v_1, v_1_m,x_1,iter_num_1,optimal_value_mark_1,optimal_value_1] = main_alin(input, -1);

 s_time_1(i) = time_1;

  s_iter_1(i) = iter_num_1;

  s_opt_1(i) = optimal_value_1(optimal_value_mark_1);
  
 s_v_1(i)=abs((s_opt_1(i) - fval)/fval);
 % s_v_1(i) = v_1(v_1_m-1);

 [time0, v0,v0_m,x0,iter_num0,optimal_value_mark0,optimal_value0] = main_alin(input, 0); 

 s_time0(i) = time0;
 
 s_iter0(i) = iter_num0;

 s_opt0(i) = optimal_value0(optimal_value_mark0);
%  
%  value_diff0 = optimal_value0(1:optimal_value_mark0) - optimal_value0(optimal_value_mark0);
s_v0(i)=abs((s_opt0(i) - fval)/fval);
 
% s_v0(i) = v0(v0_m-1);
%  hold on
   [time1,v1,v1_m,x1,iter_num1,optimal_value_mark1,optimal_value1] = main_alin(input, 1);
%   
    s_time1(i) = time1;
       s_iter1(i) = iter_num1;
%   
     s_opt1(i) = optimal_value1(optimal_value_mark1);
   s_v1(i)=abs((s_opt1(i) - fval)/fval);
%    % s_v1(i) = v1(v1_m-1);
%  hold on
%   value_diff1 = optimal_value1(1:optimal_value_mark1) - optimal_value1(optimal_value_mark1);
%s_v1(i)=log(value_diff1(optimal_value_mark1-1));
  [time2,v2, v2_m,x2,iter_num2,optimal_value_mark2,optimal_value2] = main_alin(input, 2);

 s_time2(i) = time2;

  s_iter2(i) = iter_num2;

  s_opt2(i) = optimal_value2(optimal_value_mark2);
s_v2(i)=abs((s_opt2(i) - fval)/fval);
 % s_v2(i) = v2(v2_m-1);
%value_diff2 = optimal_value2(1:optimal_value_mark2) - optimal_value2(optimal_value_mark2);
%s_v0(i)=log(value_diff2(optimal_value_mark2-1));
end
aver_runtime_1 = mean(s_time_1);
aver_runtime0 = mean(s_time0);
aver_runtime1 = mean(s_time1);
aver_runtime2 = mean(s_time2);
std_r_1 = std(s_time_1);
std_r0 = std(s_time0);
std_r1 = std(s_time1);
std_r2 = std(s_time2);
aver_iter_1 = mean(s_iter_1);
aver_iter0 = mean(s_iter0);
aver_iter1 = mean(s_iter1);
aver_iter2 = mean(s_iter2);
stdi_1 = std(s_iter_1);
stdi0 = std(s_iter0);
stdi1 = std(s_iter1);
stdi2 = std(s_iter2);
% aver_opt0 = mean(s_opt0);
% aver_opt1 = mean(s_opt1);
% aver_opt2 = mean(s_opt2);
% stdo0 = std(s_opt0);
% stdo1 = std(s_opt1);
% stdo2 = std(s_opt2);
aver_v_1 = mean(s_v_1);
aver_v0 = mean(s_v0);
aver_v1 = mean(s_v1);
aver_v2 = mean(s_v2);
stdv_1 = std(s_v_1);
stdv0 = std(s_v0);
stdv1 = std(s_v1);
stdv2 = std(s_v2);
 fprintf('selective: %.2f,%.2f,%.2f,%.2f,%.4f,%.4f\n', aver_iter2,stdi2,aver_runtime2,std_r2,aver_v2,stdv2);
  fprintf('bundle: %.2f,%.2f,%.2f,%.2f,%.4f,%.4f\n', aver_iter1,stdi1,aver_runtime1,std_r1,aver_v1,stdv1);
   fprintf('peaceman: %.2f,%.2f,%.2f,%.2f,%.4f,%.4f\n', aver_iter0,stdi0,aver_runtime0,std_r0,aver_v0,stdv0);
 fprintf('douglas: %.2f,%.2f,%.2f,%.2f,%.4f,%.4f\n', aver_iter_1,stdi_1,aver_runtime_1,std_r_1,aver_v_1,stdv_1);

   % %figure for blow up
% loglog(1:optimal_value_mark0,optimal_value0(1:optimal_value_mark0),'-^');
% %plot(1:optimal_value_mark0,optimal_value0(1:optimal_value_mark0),'-^');
% hold on;
% loglog(1:optimal_value_mark1,optimal_value1(1:optimal_value_mark1),'-*');
% %plot(1:optimal_value_mark1,optimal_value1(1:optimal_value_mark1),'-*');
% hold on;
% %plot(1:optimal_value_mark2,optimal_value2(1:optimal_value_mark2),'-o');
% loglog(1:optimal_value_mark2,optimal_value2(1:optimal_value_mark2),'-o');
% ylabel('Objective value');
% xlabel('Iterations');
% title('Change of optimal value w.r.t. the number of iterations');
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