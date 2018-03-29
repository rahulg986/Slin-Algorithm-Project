function history = rpdmm_overlap(A,b,K,rho,J,groups,gsize,overlap,alpha,lambda)

t_start = tic;

[m, n] = size(A);

MAX_ITER = 10000;

di = 2;
tK = min(K,di);

nu = 2*(J+1)./(tK.*(2*(J+1)-K));

% J = 50;
x = 0.001*eye(gsize*J,1);
z = zeros(n,1);
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

UU = -UU; % Aj^T
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

for iter = 1:MAX_ITER
    
%     idx = randperm(J,K); 
    idx_s = mod(idx_e+1,J+1);
    idx_e = idx_s+K-1;
    idx = mod(idx_s:idx_e, J+1);
     %fprintf('%d\t',idx);
     %fprintf('\n');
    pos_z = find(idx == 0);  % global variable    
    idx(pos_z) = [];
    
    % primal update
    if length(pos_z)
        zold = z;
        q = Atb + rho*(UUU.*z - Ayhat);    % temporary value
        if( m >= n )    % if skinny
           z = U \ (L \ q);
        else            % if fat
           z = q.*ud - ud.*(A'*(U \ ( L \ (A* (ud.*q) ) )));
        end
    end
    
    xold = x;
    e_idx = 0;
    for i = 1:length(idx)
        s_idx = gsize*(idx(i)-1) + 1; %s_idx{idx(i)}; %idx_e + 1; 
        e_idx = s_idx + gsize-1; %e_idx{idx(i)}; %idx_s + length(groups{i}) - 1;
        %s_idx = e_idx + 1; 
        %e_idx = s_idx + length(groups{i}) - 1;
        xyhat = x(s_idx:e_idx) - yhat(s_idx:e_idx);
        x(s_idx:e_idx) = xyhat * max( 1 - alpha(i)./norm(xyhat) , 0);
        normg = normg + rho*(norm(x(s_idx:e_idx)) - norm(xold(s_idx:e_idx)) )*alpha(idx(i));
        
        rsdl(s_idx:e_idx) = rsdl(s_idx:e_idx) + x(s_idx:e_idx) - xold(s_idx:e_idx);
    end

    if length(pos_z)    % update all dual variables
        rsdl_z = UU'*(z-zold);            % Aj*xj
        rsdl = rsdl + rsdl_z;
    end
    yhat = nu*rsdl;
    
    Ayhat = UU*yhat;
    
    % obj
    history.objval(iter) = (norm(A*z - b,2)^2/(2*gamma) + normg);
    xzold = xz;
    xz = [x(:);z];
    
    history.r_norm(iter)= norm(xzold - xz)/norm(xzold);
    history.s_norm(iter)= norm(rsdl);
    
    history.time(iter) = toc(t_start);
    
      fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n', iter, ...
              history.r_norm(iter),history.s_norm(iter), history.objval(iter),history.time(iter));
    
    if history.r_norm(iter)+ history.s_norm(iter) < 1e-4
        break;
    end
    
end % end iter

%history.x = x;
%history.z = z;

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