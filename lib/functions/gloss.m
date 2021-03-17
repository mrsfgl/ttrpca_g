function [L, S, Nt, obj_val] = gloss(Y, param)
% [L, S, Nt, obj_val] = gloss(Y, param)
% Graph Regularized Low rank plus Smooth-Sparse Decomposition

N = ndims(Y);
sz = size(Y);
max_iter = param.max_iter;
err_tol = param.err_tol;
alpha = param.alpha;
theta = param.theta;
lambda = param.lambda;
gamma = param.gamma;
psi = param.psi;
beta_1 = param.beta_1;
beta_2 = param.beta_2;
beta_3 = param.beta_3;
beta_4 = param.beta_4;
L = cell(1, N);
for i=1:N
    L{i} = zeros(size(Y));
end
La = L;
S = zeros(size(Y));
Nt = zeros(size(Y));
W = zeros(size(Y));
Z = zeros(size(Y));
Phi = get_graphL(Y, 5);
inv_Phi = ((theta/beta_4)*Phi+eye(size(Phi)))^-1;
D = convmtx([1,-1], size(Y,1));
D(:,end) = [];
D(end,1) = -1;
if beta_3~=0 && beta_2~=0
    invD = (beta_3*eye(sz(1))+beta_2*(D'*D))^-1;
else
    invD = zeros(sz(1));
end
Lam{1} = cell(1, N);
for i=1:N
    Lam{1}{i} = zeros(size(Y));
end
Lam{4} = Lam{1};
Lam{2} = zeros(size(Y));
Lam{3} = zeros(size(Y));

timeL = [];
timeS = []; 
timeN = [];
timeZ = [];
timeW = [];
timeDual = [];
iter = 1;
obj_val = compute_obj(Y,L,La,S,Nt,W,Z,Lam,D,Phi,param);
while true
    %% L Update
    tstart = tic;
    [L, ~] = soft_hosvd_wgr(Y-S-Nt, Lam{1}, La, Lam{4}, psi, beta_1, beta_4);
    timeL(end+1)=toc(tstart);
    %% La Update
    La = graph_reg_update(L,Lam{4},inv_Phi);
    %% S Update
    tstart = tic;
    temp1 = zeros(size(Y));
    for i=1:N
        temp1 = temp1+beta_1*(Y-L{i}-Nt+Lam{1}{i});
    end
    temp2 = beta_3*(W+Lam{3});
    Sold = S;
    S = soft_threshold((temp1+temp2), lambda)./(N*beta_1+beta_3);
    timeS(end+1)=toc(tstart);
    %% N update
    tstart = tic;
    Nt = zeros(size(Y));
    for i=1:N
        Nt = Nt+(beta_1/(N*beta_1+alpha)).*(Y+Lam{1}{i}-L{i}-S);
    end
    timeN(end+1)=toc(tstart);
    %% W Update
    tstart = tic;
    W = invD*(beta_3*Runfold(S-Lam{3})+beta_2*D'*Runfold(Z+Lam{2}));
    W = reshape(W, sz);
    timeW(end+1)=toc(tstart);
    
    %% Z Update
    tstart = tic;
    Z = soft_threshold(mergeTensors(D, W, 2, 1)-Lam{2}, gamma/(beta_2+eps));
    timeZ(end+1)=toc(tstart);
    %% Dual Updates
    tstart = tic;
    temp = 0;
    for i=1:N
        Lam1_up = Y-L{i}-S-Nt;
        temp = temp + norm(Lam1_up(:))^2;
        Lam{1}{i} = Lam{1}{i}+Lam1_up;
        Lam{4}{i} = Lam{4}{i}-(L{i}-La{i});
    end
    temp = sqrt(temp)/(sqrt(N)*norm(Y(:)));
    Lam{2} = Lam{2}-mergeTensors(D, W, 2, 1)+Z;
    Lam{3} = Lam{3}-S+W;
    timeDual(end+1)=toc(tstart);
    
    %% Error and objective calculations
    obj_val(iter+1) = compute_obj(Y,L,La,S,Nt,W,Z,Lam,D,Phi,param);
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

function [val, term] = compute_obj(Y,L,La,S,Nt,W,Z,Lam,D,Phi,param)

alpha = param.alpha;
theta = param.theta;
lambda = param.lambda;
gamma = param.gamma;
beta_1 = param.beta_1;
beta_2 = param.beta_2;
beta_3 = param.beta_3;
beta_4 = param.beta_4;
N = length(L);
term = zeros(1,9);
for i=1:N
    term(6) = term(6) + beta_1/2*sum((Y-S-L{i}-Nt+Lam{1}{i}).^2,'all');
    term(1) = term(1) + comp_nuclear(L{i}, i);
    term(9) = term(9) + beta_4/2*sum((L{i}-La{i}-Lam{4}{i}).^2,'all');
    term(2) = term(2) + theta*comp_gr_reg(La{i}, Phi);
end
term(7) = (beta_2/2)*sum((mergeTensors(D, W, 2, 1)-Lam{2}-Z).^2,'all');
term(8) = (beta_3/2)*sum((S-W-Lam{3}).^2,'all');
term(3) = lambda*sum(abs(S),'all');
term(5) = (alpha/2)*norm(Nt(:))^2;
term(4) = gamma*sum(abs(Z),'all');
val = sum(term);
end