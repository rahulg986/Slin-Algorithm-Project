function history = sadmm_overlap(A,b,K,rho,J,groups,alpha,lambda)

t_start = tic;

[m, n] = size(A);

MAX_ITER = 1000;

di = J+1;
tK = min(K,di);

tau = 1*K./(tK.*(2*(J+1)-K));
nu = 1./tK;

% J = 50;
x = 0.001*eye(100,J);
z = zeros(n,1);
y = zeros(100*J,1);
yhat = zeros(100*J,1);

% save a matrix-vector multiply
gamma = lambda*J;
Atb = A'*b/gamma;

UU = speye(n,100);
for i = 1:length(groups)-1
   UU = [UU,circshift(eye(n,100),90*i)];
end

UUU = full(sum(UU,2));

UU = -UU; % Aj^T
ud = 1./(UUU*rho);

[L U] = factor(A, gamma, ud);
%save LU.mat L U;
%load LU.mat;

alpha = alpha/rho;
xz = 0;
Ayhat = 0;

for iter = 1:MAX_ITER
        
    % primal update
    zold = z;
    q = Atb + rho*(UUU.*z - Ayhat);    % temporary value
    if( m >= n )    % if skinny
       z = U \ (L \ q);
    else            % if fat
       z = q.*ud - ud.*(A'*(U \ ( L \ (A* (ud.*q) ) )));
    end
    
    xold = x;
    xyhat = x - reshape(yhat,[100,J]);
    x = xyhat * diag( max( 1 - alpha'./sqrt(sum(xyhat.^2,1)) , 0) );
    normg = rho*sqrt(sum(x.^2,1))*alpha;
    
    % dual update 
    % no global variable, just update its own dual variable
    yold = y;
    rsdl = x(:) + UU'*z;
    y = y + tau*rsdl;
    yhat = y + nu*rsdl;
    
    Ayhat = UU*yhat;
    
    % obj
    history.objval(iter) = norm(A*z - b,2)^2/(2*gamma) + normg;
    xzold = xz;
    xz = [x(:);z];
    
    history.r_norm(iter)= norm(xzold - xz)/norm(xzold);
    history.s_norm(iter)= norm(yold - y)/norm(yold);
    
    history.time(iter) = toc(t_start);
    
%     fprintf('%3d\t%10.4f\t%10.4f\t%10.2f\t%10.4f\n', iter, ...
%             history.r_norm(iter),history.s_norm(iter), history.objval(iter),history.time(iter));
    
    if history.r_norm(iter)+ history.s_norm(iter) < 1e-4
        break;
    end
    
end % end iter

end % end pdmm func

function [L U] = factor(A, gamma, ud)
    [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( A'*A/gamma + diag(1./ud), 'lower' );
    else            % if fat
       L = chol( gamma*speye(m) + A*diag(ud)*A', 'lower' );
    end

    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
end