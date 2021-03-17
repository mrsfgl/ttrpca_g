function [L, S, obj_val] = horpca(Y, param)
% [L, S, Nt, obj_val] = horpca(Y, param)
% Higher order Robust PCA

N = ndims(Y);

max_iter = param.max_iter;
err_tol = param.err_tol;
lambda = param.lambda;
psi = param.psi;

beta_1 = param.beta_1;

L = cell(1, N);
for i=1:N
    L{i} = zeros(size(Y));
end
S = zeros(size(Y));

Lam = cell(1, N);
for i=1:N
    Lam{i} = zeros(size(Y));
end

timeL = [];
timeS = []; 
timeDual = [];
iter = 1;
obj_val = compute_obj(Y,L,S,Lam,param);
while true
    %% L Update
    tstart = tic;
    [L, ~] = soft_hosvd(Y-S, Lam, psi, 1/beta_1);
    timeL(end+1)=toc(tstart);
    %% S Update
    tstart = tic;
    temp1 = zeros(size(Y));
    for i=1:N
        temp1 = temp1+(L{i}-Lam{i});
    end
    Sold = S;
    S = soft_threshold(beta_1*(Y-temp1/N), lambda)./(beta_1);
    timeS(end+1)=toc(tstart);
    %% Dual Updates
    tstart = tic;
    temp = 0;
    for i=1:N
        Lam1_up = Y-L{i}-S;
        temp = temp + norm(Lam1_up(:))^2;
        Lam{i} = Lam{i}+Lam1_up;
    end
    temp = sqrt(temp)/(sqrt(N)*norm(Y(:)));
    timeDual(end+1)=toc(tstart);
    
    %% Error and objective calculations
    obj_val(iter+1) = compute_obj(Y,L,S,Lam,param);
%     err = abs(obj_val(iter)-obj_val(iter+1))/obj_val(iter);
    err = max(norm(S(:)-Sold(:))/norm(Sold(:)), temp);
    iter = iter+1;
    
    if err<=err_tol
        disp('Converged!')
        break;
    end
    if iter>=max_iter
        disp('Max iter')
        break;
    end
end
temp = zeros(size(Y));
for i=1:N
    temp = temp+L{i};
end
L = temp/N;

end

function [val, term] = compute_obj(Y,L,S,Lam,param)

lambda = param.lambda;
beta_1 = param.beta_1;
N = length(L);
term = zeros(1,3);
for i=1:N
    term(3) = term(3) + beta_1/2*sum((Y-S-L{i}+Lam{i}).^2,'all');
    term(1) = term(1) + comp_nuclear(L{i}, i);
end
term(2) = lambda*sum(abs(S),'all');
val = sum(term);
end