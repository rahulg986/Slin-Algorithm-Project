function history = rbsumm_overlap(A,b,rho,J,groups,gsize,overlap,alpha,lambda)

t_start = tic;

[m, n] = size(A);

MAX_ITER = 10000;


% J = 50;
x = 0.001*eye(gsize*J,1);
z = zeros(n,1);
y = zeros(gsize*J,1);
yhat = zeros(gsize*J,1);

% save a matrix-vector multiply
gamma = lambda*J;
Atb = A'*b/gamma;

%UU = zeros(n);

UU = speye(n,gsize);
for i = 1:length(groups)-1
   UU = [UU,circshift(eye(n,gsize),(gsize-overlap)*i)];
end

UUU = full(sum(UU,2));

ud = 1./(UUU*rho);

[L U] = factor(A, gamma, ud);
%save LU.mat L U;
%load LU.mat;

alpha = alpha/rho;
normg = rho*norm(x)*alpha(1);
rsdl = x;
xz = 1;
Ayhat = 0;
idx_e = -1;
yold = 0.000001;
uz = UU'*z;

for iter = 1:MAX_ITER
    
%     idx = randperm(J,K); 
    idx = randi(J+2)-1;
    %idx = mod(iter-1,J+2);
    
    switch idx
        case 0 
    % primal update
            zold = z;
            q = Atb + rho*(UUU.*z - Ayhat);    % temporary value
            if( m >= n )    % if skinny
               z = U \ (L \ q);
            else            % if fat
               z = q.*ud - ud.*(A'*(U \ ( L \ (A* (ud.*q) ) )));
            end
            uz = UU'*z;
        case J+1
            % dual update 
            yold = y;
            ii = 2;
            eta = ii^2+1/(ii+iter);
            y = y + eta*(x+uz);
            Ayhat = UU*y;
        otherwise        
    
            idx_s = gsize*(idx-1) + 1;
            idx_e = idx_s + length(groups{idx}) - 1;
            
            xold = x;
            yz = uz(idx_s:idx_e) - y(idx_s:idx_e);
      
            x(idx_s:idx_e) = yz * max( 1 - alpha(i)/norm(yz) , 0);
%            normg = normg + rho*( norm(x(idx_s:idx_e)) - norm(xold(idx_s:idx_e)) )*alpha(i);

    end
    
    % obj
    normg = 0;
    idx_e = 0;
    for i = 1:length(groups)
        idx_s = idx_e + 1; 
        idx_e = idx_s + length(groups{i}) - 1;
        normg = normg + rho*sqrt(sum(x(idx_s:idx_e).^2,1))*alpha(i);
    end
    history.objval(iter) = (norm(A*z - b,2)^2/(2*gamma) + normg);
    xzold = xz;
    xz = [x(:);z];
    
    history.r_norm(iter)= norm(xzold - xz)/norm(xzold);
    history.s_norm(iter)= norm(yold - y)/norm(yold);
    
    history.time(iter) = toc(t_start);
    
    fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n', iter, ...
            history.r_norm(iter),history.s_norm(iter), history.objval(iter),history.time(iter));
    
    if history.r_norm(iter)+ history.s_norm(iter) < 1e-4
        break;
    end
    
end % end iter

history.x = x;
history.z = z;

end % end pdmm func

function [L U] = factor(A, gamma, ud)
    [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( A'*A/gamma + diag(1./ud), 'lower' );
    else            % if fat
       T = gamma*speye(m) + (A.*repmat(ud',m,1))*A';
       L = chol( T, 'lower' );
    end

    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
end