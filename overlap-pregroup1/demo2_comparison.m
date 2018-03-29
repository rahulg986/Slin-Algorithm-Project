%clear;
%clc;

%% Prepare the data
ng=50;  % number of groups 
g_size=100;  % group size
overlap=10;  % number of variables overlapped between two consecutive groups
n=1000;      % sample size
p=ng*(g_size-overlap)+overlap;   % totoal number of variables 
input.EPI = 1e-5;
input.n = n;
input.p = p;
input.ng = ng;
input.g_size = g_size;
input.overlap = overlap;

%generate toy data
% A: design matrix
% b: output
% T: ng \times p, indicate each group contains which variables 
% Tx: weight for each group
% x: true regression coefficients

 r_T=kron(1:ng, ones(1,g_size));
    c_T=zeros(ng*g_size, 1);
    s=1;
    for g=1:ng
        c_T((g-1)*g_size+1:g*g_size,1)=s:s+g_size-1;
        s=s+g_size-overlap;
    end
    T=sparse(r_T, c_T, 1, ng, p);
     Tx=ones(ng,1);  % uniform weight
     tmp=1:p;
    x=((-1).^tmp(:)).*(exp(-(tmp(:)-1)/100));  
    sn2 = 1;  % signal to noise ratio
    
  %  A =randn(n,p);
   % b= A*x+sn2*rand(n,1);     
input.A=A;
input.b=b;

input.T = T;
%tau=0.1*norm(A'*b,Inf);
lamda1 = 1/g_size;
input.lambda1 = lamda1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Call the three methods for comparison
% [x0,iter_num0,optimal_value_mark0,optimal_value0] = main_alin(input, 0); 
% hold onb_h(2)= b_h(2) + yhi(i,z_h(:,2),lambda1,T) - a_h(:,2)'*z_h(:,2);
%[x1,iter_num1,optimal_value_mark1,optimal_value1] = main_alin(input, 1);
%  hold on
% value_diff1 = optimal_value1(1:optimal_value_mark1) - optimal_value1(optimal_value_mark1);
[x2,iter_num2,optimal_value_mark2,optimal_value2] = main_alin(input, 2);
%value_diff2 = optimal_value2(1:optimal_value_mark2) - optimal_value2(optimal_value_mark2);

% figure for blow up
% plot(1:optimal_value_mark0,optimal_value0(1:optimal_value_mark0),'-^');
% hold on;
% %plot(1:optimal_value_mark1,optimal_value1(1:optimal_value_mark1),'-*');
% %hold on;
% plot(1:optimal_value_mark2,optimal_value2(1:optimal_value_mark2),'-o');
% ylabel('optimal value');
% xlabel('Iterations');
% title('Change of optimal value w.r.t. the number of iterations');
%%%%%%%%figure ln(f(x)-f(x*))
 %figure
%value_diff = optimal_value2(1:optimal_value_mark2) - optimal_value2(optimal_value_mark2);
% %if (update_selector ==  1)
%     plot(1:(optimal_value_mark1-1), log(value_diff1(1:optimal_value_mark1-1)), '-*');
%     hold on;
% %else
  %   plot(1:(optimal_value_mark2-1), log(value_diff2(1:optimal_value_mark2-1)), '-o');
% %end
% plot(1:(optimal_value_mark-1), log(value_diff(1:optimal_value_mark-1)), '-*');
%  xlabel('Iterations');
%  ylabel('ln(F(x) - F(x*))');
%  title('Change of ln(F(x) - F(x^*)) w.r.t. the number of iterations');


%% Before figure, let's output the value first to get some sense
%scatter(1:10,iter_plot);
%xlabel('-log_{10}(eps)')
%ylabel('Total number of steps')
%title_str = strcat('m=', num2str(n), ', n=', num2str(p), ', eps=', num2str(input.EPI));
%title(title_str)