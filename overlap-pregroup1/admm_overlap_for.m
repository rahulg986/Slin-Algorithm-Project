function history = admm_overlap_for(A,b,K,rho,J,groups,gsize,overlap,alpha,lambda)

t_start = tic;

[m, n] = size(A);

MAX_ITER = 1000;

% J = 50;
%x = zeros(m,1);
z = randn(n,1);
%y = zeros(m,1);

% save a matrix-vector multiply
gamma = lambda*J;
Atb = A'*b/gamma;

UU{1} = speye(n,gsize);
UUU = zeros(n,1);
UUU(groups{1}) = 1;
x{1} = zeros(length(groups{1}),1);
y{1} = zeros(length(groups{1}),1);
for i = 2:length(groups)
   UU{i} = circshift(UU{i-1},gsize-overlap);
   %UUU = UUU + [UU{i}, zeros(n,)];
   UUU(groups{i}) = UUU(groups{i}) + 1;
   x{i} = zeros(length(groups{i}),1);
   y{i} = zeros(length(groups{i}),1);
end

% UU = speye(n,100);
% for i = 1:length(groups)-1
%    UU = [UU,circshift(eye(n,100),90*i)];
% end

% UUU = full(sum(UU,2));

ud = 1./(UUU*rho);

[L U] = factor(A, gamma, ud);
%save LU.mat L U;
%load LU.mat;

alpha = alpha/rho;

for iter = 1:MAX_ITER
    
    % primal update
    xold = x;
    yold = y;
    normg = 0;
    uxy = 0;
    history.s_norm(iter) = 0;
    for i = 1:length(groups)
        y{i} = y{i} + x{i} - UU{i}'*z;
        
        yz = UU{i}'*z - y{i};
        x{i} = yz * diag( max( 1 - alpha(i)/norm(yz) , 0) );
        normg = normg + rho*norm(x{i})*alpha(i);
        
        uxy = uxy + UU{i}*(x{i} + y{i});
        
        history.s_norm(iter)= history.s_norm(iter) + norm(yold{i} - y{i})/norm(yold{i});
    end
    
    zold = z;
    q = Atb + rho*uxy; 
    if( m >= n )    % if skinny
       z = U \ (L \ q);
    else            % if fat
       z = q.*ud - ud.*(A'*(U \ ( L \ (A* (ud.*q) ) )));
    end
    
    % dual update 
%     yold = y;
%     y = y + x - UU'*z;
    
    % obj
    history.objval(iter) = norm(A*z - b,2)^2/(2*gamma) + normg;
    %xzold = xz;
    %xz = [x(:);z];
    
    history.r_norm(iter)= norm(zold - z)/norm(zold);
    %history.s_norm(iter)= norm(yold - y)/norm(yold);
    
    history.time(iter) = toc(t_start);
    
    fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n', iter, ...
            history.r_norm(iter),history.s_norm(iter), history.objval(iter),history.time(iter));
    
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