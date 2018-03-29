function history = rpdmm_overlap(A,b,K,rho,J,groups,alpha,lambda)

t_start = tic;

[m, n] = size(A);

MAX_ITER = 5000;

di = 2;
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

%UU = zeros(n);

UU = speye(n,100);
for i = 1:length(groups)-1
   UU = [UU,circshift(eye(n,100),90*i)];
end

UUU = full(sum(UU,2));

UU = -UU; % Aj^T
ud = 1./(UUU*rho);

%[L U] = factor(A, gamma, ud);
%save LU.mat L U;
load LU.mat;

alpha = alpha/rho;
normg = rho*sqrt(sum(x.^2,1))*alpha;
rsdl = x;
xz = 0;
Ayhat = 0;
idx_e = -1;

for iter = 1:MAX_ITER
    
%    idx = randperm(J+1,K);
%     pos_z = find(idx == J+1); 
    
     idx_s = mod(idx_e+1,J+1);
     idx_e = idx_s+K-1;
     idx = mod(idx_s:idx_e, J+1);
    pos_z = find(idx == 0);  % global variable
    
    idx(pos_z) = [];
    %pos_x = find(idx ~= J+1);
    
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
    
    if length(idx)
        xtold = x(:,idx);
        yhat_t = reshape(yhat,[100,J]);
        xyhat = xtold - yhat_t(:,idx);
        xt = xyhat * diag( max( 1 - alpha(idx)'./sqrt(sum(xyhat.^2,1)) , 0) );
        x(:,idx) = xt;
        normg = normg + rho*( sqrt(sum(xt.^2,1)) - sqrt(sum(xtold.^2,1)) )*alpha(idx);
        
        rsdl(:,idx) = rsdl(:,idx) + xt - xtold;
    end
    %     normg = 0;
    %     xold = x;
    %     idx_e = 0;
    %     for i = 1:length(idx)
    %         idx_s = idx_e + 1; 
    %         idx_e = idx_s + length(groups{i}) - 1;
    %         xyhat = x(idx_s:idx_e) - yhat(idx_s:idx_e);
    %         x(idx_s:idx_e) = xyhat * diag( max( 1 - alpha(i)./sqrt(sum(xyhat.^2,1)) , 0) );
    %         normg = normg + rho*sqrt(sum(x(idx_s:idx_e).^2,1))*alpha(i);
    %     end

        % dual update 
        % no global variable, just update its own dual variable
    yold = y;
    if length(pos_z)    % update all dual variables
        rsdl_z = UU'*(z-zold);            % Aj*xj
        rsdl = rsdl + reshape(rsdl_z,[100,J]);
    end
    y = y + tau*rsdl(:);
    yhat = y + nu*rsdl(:);
    
%     yold = y;
%     rsdl = x(:) + UU'*z;
%     y = y + tau*rsdl;
%     yhat = y + nu*rsdl;
    
    Ayhat = UU*yhat;
    
    % obj
    history.objval(iter) = norm(A*z - b,2)^2/(2*gamma) + normg;
    xzold = xz;
    xz = [x(:);z];
    
    history.r_norm(iter)= norm(xzold - xz)/norm(xzold);
    history.s_norm(iter)= norm(yold - y)/norm(yold);
    
    history.time(iter) = toc(t_start);
    
    fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n', iter, ...
            history.r_norm(iter),history.s_norm(iter), history.objval(iter),history.time(iter));
    
    if history.r_norm(iter)+ history.s_norm(iter) < 1e-3
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
       L = chol( gamma*speye(m) + A*diag(ud)*A', 'lower' );
    end

    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
end