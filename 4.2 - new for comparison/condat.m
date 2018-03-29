% 
% The iterative over-relaxed Chambolle-Pock algorithm converges to 
% the image x minimizing 1/2*||Ax-b||^2/2 + lambda1*|x|_1 +
% lambda2*|x_{j+1} - x_j|
% There is only one parameter to tune, tau, which controls the convergence speed.
%
%Code written by Laurent Condat, CNRS research fellow in 
% Dept. of Images and Signals of GIPSA-lab, Univ. Grenoble Alpes, Grenoble, France.
%?
% See the description of the algorithm in 
% L. Condat, "A primal-dual splitting method for convex optimization 
% involving Lipschitzian, proximable and linear composite terms", 
% J. Optimization Theory and Applications, vol. 158, no. 2, pp. 460-479, 
% 2013



function main
tic
	Nbiter= 300;
	lambda = 25; 
   % lambda2 = 25;
	tau = 0.05;
	eta = 0.99;
%     n=300;
% p=1000;
% input.EPI = 1e-3;
% input.n = n;
% input.p = p;
% x=zeros(p,1);
% pos=unidrnd(p,160,1);
% for i=1:numel(pos)
%     if (mod(pos(i),2)==1) 
%         x(pos(i))=1;
%     else
%         x(pos(i))=-1;
%     end
% end
% R = zeros(p-1,p);
% for i = 1:(p-1)
%     R(i,i) = -1;
%     R(i,i+1) = 1;
% end
% 
% input.A=rand(n,p);
% A = input.A;
% input.b=A*x + normrnd(0,10^(-4),size(A,1),1);
% b = input.b;
% tau=0.1*norm(A'*b,Inf);
% input.lambda1 = tau;
% input.lambda2 = tau;
	
	y  = double(imread('parrotgray.png'));   % Initial image
	figure(1);
	imshow(y/255);
	rng(0);
	y = y+randn(size(y))*30; % White Gaussian noise is added to the image
	figure(2);
	imshow(y/255);
	
	opD = @(x) cat(3,[diff(x,1,1);zeros(1,size(x,2))],[diff(x,1,2) zeros(size(x,1),1)]);
	opDadj = @(u) -[u(1,:,1);diff(u(:,:,1),1,1)]-[u(:,1,2) diff(u(:,:,2),1,2)];	
	prox_tau_f = @(x) (x+tau*y)/(1+tau);
	prox_sigma_g_conj = @(u) bsxfun(@rdivide,u,max(sqrt(sum(u.^2,3))/lambda,1));
	
	x2 = y; 		% Initialization of the solution
	u2 = prox_sigma_g_conj(opD(x2));	% Initialization of the dual solution
	sigma = 1/tau/8;
%a = x_k - tau*
		
	for iter = 1:Nbiter
		x = prox_tau_f(x2-tau*opDadj(u2));
		u = prox_sigma_g_conj(u2+sigma*opD(2*x-x2));
		x2 = x+eta*(x-x2);
		u2 = u+eta*(u-u2);
		if mod(iter,20)==0
			fprintf('nb iter:%4d  cost value: %f\n',iter,...
				norm(x-y,'fro')^2/2+lambda*sum(sum(sqrt(sum(opD(x).^2,3)))));
			figure(3);
			imshow(x/255);
			colormap gray
		end
	end
	figure(3);
	imshow(x/255);
toc	
end
