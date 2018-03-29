%conjugate gradient method with no precondition
function [ x_cur,i ] =prec_conjug( A,f,rho,x_cur )
i=0;
dev_cur=(A'*(A*x_cur)+rho.*x_cur)-f;
d_cur=-dev_cur;
norm_d=norm(d_cur);
norm_dev=norm_d;
eps=10^(-5);
%while (delta_new > eps^2*delta_0)
while (norm_dev > eps)
   %%%%%q is basicly Q*d_cur
    q=A'*(A*d_cur) + rho.*d_cur;
    %%%%Step length and new points
    tau=norm_dev^2/(d_cur'*q);
    x_new=x_cur + tau*d_cur;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    dev_diff=tau*q;
    dev_new=dev_cur + dev_diff;
    alpha=dev_new'*dev_diff/norm_dev^2;
    d_new=-dev_new + alpha*d_cur;
    d_cur=d_new;
    norm_dev=norm(dev_new);
    x_cur=x_new;
    dev_cur=dev_new;
    i=i+1;
end


end

