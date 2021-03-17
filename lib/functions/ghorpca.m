function [L, S, N, obj_val] = ghorpca(Y, param)
% [L, S, Nt, obj_val] = gloss(Y, param)
% Graph Regularized Low rank plus Smooth-Sparse Decomposition

N = ndims(Y);
sz = size(Y);
max_iter = param.max_iter;
err_tol = param.err_tol;
alpha = param.alpha;
theta = param.theta;
lambda = param.lambda;
psi = param.psi;
beta_1 = param.beta_1;
beta_2 = param.beta_4;
L = cell(1, N);
for i=1:N
    L{i} = zeros(size(Y));
end
La = L;
S = zeros(size(Y));
Nt = zeros(size(Y));
Phi = get_graphL(Y, 3, true);
inv_Phi = ((theta/beta_2)*Phi+eye(size(Phi)))^-1;
Lam{1} = cell(1, N);
for i=1:N
    Lam{1}{i} = zeros(size(Y));
end
Lam{2} = Lam{1};

timeL = [];
timeS = []; 
timeN = [];
timeDual = [];
iter = 1;
obj_val = compute_obj(Y,L,La,S,Nt,Lam,Phi,param);
while true
    %% L Update
    tstart = tic;
    [L, ~] = soft_hosvd_wgr(Y-S-Nt, Lam{1}, La, Lam{2}, psi, beta_1, beta_2);
    timeL(end+1)=toc(tstart);
    %% La Update
    La = graph_reg_update(L,Lam{2},inv_Phi);
    %% S Update
    tstart = tic;
    temp1 = zeros(size(Y));
    for i=1:N
        temp1 = temp1+L{i}-Lam{1}{i};
    end
    Sold = S;
    S = soft_threshold(beta_1*(N*(Y-Nt)-temp1), lambda)./(N*beta_1);
    timeS(end+1)=toc(tstart);
    %% N update
    tstart = tic;
    Nt = zeros(size(Y));
    for i=1:N
        Nt = (beta_1/(N*beta_1+alpha)).*(Y+Lam{1}{i}-L{i}-S);
    end
    timeN(end+1)=toc(tstart);
    %% Dual Updates
    tstart = tic;
    temp = 0;
    for i=1:N
        Lam1_up = Y-L{i}-S-Nt;
        temp = temp + norm(Lam1_up(:))^2;
        Lam{1}{i} = Lam{1}{i}+Lam1_up;
        Lam{2}{i} = Lam{2}{i}-(L{i}-La{i});
    end
    temp = sqrt(temp)/(sqrt(N)*norm(Y(:)));
    timeDual(end+1)=toc(tstart);
    
    %% Error and objective calculations
    obj_val(iter+1) = compute_obj(Y,L,La,S,Nt,Lam,Phi,param);
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

function [val, term] = compute_obj(Y,L,La,S,Nt,Lam,Phi,param)

alpha = param.alpha;
theta = param.theta;
lambda = param.lambda;
beta_1 = param.beta_1;
beta_2 = param.beta_4;
N = length(L);
term = zeros(1,6);
for i=1:N
    term(5) = term(5) + beta_1/2*sum((Y-S-L{i}-Nt+Lam{1}{i}).^2,'all');
    term(1) = term(1) + comp_nuclear(L{i}, i);
    term(6) = term(6) + beta_2/2*sum((L{i}-La{i}-Lam{2}{i}).^2,'all');
    term(2) = term(2) + theta*comp_gr_reg(La{i}, Phi);
end
term(3) = lambda*sum(abs(S),'all');
term(4) = (alpha/2)*norm(Nt(:))^2;
val = sum(term);
end