%n is number of oberservations
% p is number of attributes
% clc;
% clear;
n=2^10;
p=2^14;
beta=zeros(p,1);
%value=unifrnd(-5,5,10);
pos=unidrnd(p,160,1);
for i=1:numel(pos)
    if (mod(pos(i),2)==1) beta(pos(i))=1;
    else
        beta(pos(i))=-1;
    end
end
%Generate data from multivariate normal random distribution
X=rand(n,p);
%X=[ones(n+3*n,1),X];
%X = bsxfun(@rdivide, X, sqrt(sum(X.^2,1))); % normalize predictors
% new_X=X(n+1:end,:);
% X(n+1:end,:)=[];
%data normalize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:size(X,2)
% X(:,i)=(X(:,i)-mean(X(:,i)))/std(X(:,i));
% end
% for i=1:size(new_X,2)
% new_X(:,i)=(new_X(:,i)-mean(new_X(:,i)))/std(new_X(:,i));
% end
%Generate response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=X*beta + normrnd(0,10^(-4),size(X,1),1);
%new_y=new_X*beta + normrnd(0,1,size(new_X,1),1);
tau=0.1*norm(X'*y,Inf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save('data_500.mat','X','y','new_X','new_y','beta','pos','value')
