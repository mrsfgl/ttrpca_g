function [L, S, times, obj_val, terms] = logss(Y, param)
% [L, S] = logss(Y, param)
% Low rank On Graphs plus Smooth-Sparse Decomposition

N = ndims(Y);
sz = size(Y);
mask_Y = ones(sz)>0;
mask_Y(param.ind_m) = ~mask_Y(param.ind_m);
comp_obj = param.disp;
max_iter = param.max_iter;
err_tol = param.err_tol;
lambda = param.lambda;
gamma = param.gamma;
theta = param.theta;
beta_1 = param.beta_1;
beta_2 = param.beta_2;
beta_3 = param.beta_3;
beta_4 = param.beta_4;
L = zeros(sz);
S = zeros(sz);
W = zeros(sz);
Z = zeros(sz);
U = cell(1, N);
inv_GR = cell(1, N);
Sig = U;
sz_lo = U;
p.paramn = struct;
for i = 1:N
    p.paramn.k = min(10, sz(i)-1);
    Gr = gsp_nn_graph(t2m(Y,i), p.paramn);
    [V, E] = eig(full(Gr.L));
    a = diag(E);
    ind = find(a(1:end-1)./a(2:end)>.9);
    ind = max(ind(1),round(sz(i)/2));
    U{i} = V(:,1:ind);
    Sig{i} = diag(diag(E(1:ind,1:ind)));
    inv_GR{i} = (2*theta*Sig{i}/beta_4+eye(ind))^-1;
    sz_lo{i} = [sz(1:i-1), ind, sz(i+1:N)];
end
Gx = cell(1, N);
for i=1:N
    Gx{i} = zeros(sz);
end
D = convmtx([1,-1], size(Y,1));
D(:,end) = [];
D(end,1) = -1;
if beta_3~=0 && beta_2~=0
    invD = (beta_3*eye(sz(1))+beta_2*(D'*D))^-1;
else
    invD = zeros(sz(1));
end
Lam{4} = cell(1, N);
Gp = cell(1,N);
for i=1:N
    Lam{4}{i} = zeros(sz);
    Gp{i} = zeros(sz);
end
Lam{1} = zeros(sz);
Lam{2} = zeros(sz);
Lam{3} = zeros(sz);

times = [];
iter = 1;
while true
    %% L Update
    tstart = tic;
    T1 = Y-S+Lam{1};
    T2 = zeros(sz);
    for n=1:N
        T2 = T2 + Gp{n} + Lam{4}{n};
    end
    L(mask_Y) = (beta_1*T1(mask_Y)+beta_4*T2(mask_Y))/(beta_1+N*beta_4);
    L(~mask_Y) = T2(~mask_Y)/N;
    times(iter,1) = toc(tstart);
        
    %% G update
    tstart = tic;
    for i=1:N
        Gx{i} = m2t(inv_GR{i}*(U{i}'*t2m((L-Lam{4}{i}),i)), sz_lo{i}, i);
    end
    times(iter,2) = toc(tstart);
    
    %% S Update
    tstart = tic;
    T1 = beta_1*(Y-L+Lam{1});
    T2 = W+Lam{3};
    Sold = S;
    S(mask_Y) = soft_threshold(T1(mask_Y)+beta_3.*T2(mask_Y), lambda)./(beta_1+beta_3);
    S(~mask_Y) = soft_threshold(T2(~mask_Y), lambda)./(beta_3);
    times(iter,3) = toc(tstart);
    
    %% W Update
    tstart = tic;
    W = invD*(beta_3*Runfold(S-Lam{3})+beta_2*D'*Runfold(Z+Lam{2}));
    W = reshape(W, sz);
    times(iter,4) = toc(tstart);
    
    %% Z Update
    tstart = tic;
    Z = soft_threshold(mergeTensors(D, W, 2, 1)-Lam{2}, gamma/(beta_2+eps));
    times(iter,5) = toc(tstart);
    
    %% Dual Updates
    tstart = tic;
    temp = 0;
    Lam1_up = Y-L-S;
    Lam{1}(mask_Y) = Lam{1}(mask_Y)+Lam1_up(mask_Y);
    Lam{2} = Lam{2}-mergeTensors(D, W, 2, 1)+Z;
    Lam{3} = Lam{3}-S+W;
    for n=1:N
        Gp{n} = tmprod(Gx{n}, U{n}, n);
        Lam{4}{n} = Lam{4}{n} - (L - Gp{n});
    end
    temp = sqrt(temp)/(norm(Y(:)));
    times(iter,6) = toc(tstart);
    
    %% Error calculations
    err = max(norm(S(:)-Sold(:))/norm(Sold(:)), temp);
    iter = iter+1;
    
    if comp_obj
        [obj_val(iter), terms(:,iter)] = obj_fun(Y, L, Gx, Gp, S, W, Z, Lam, Sig, mask_Y, D, param);
    end
    if err<=err_tol
        disp('Converged!')
        break;
    end
    if iter>max_iter
        disp('Max iter')
        break;
    end
end

end

function [obj_val, terms] = obj_fun(Y, L, Gx, Gp, S, W, Z, Lam, Sig, mask, D, param)

N = ndims(L);

terms(7) = 0;
for n=1:N
    temp = ndim_unfold(Gx{n},n);
    terms(1) = terms(1)+param.theta*trace(temp*temp'*Sig{n});
    terms(7) = terms(7)+param.beta_4*norm(T2V(L-Gp{n}-Lam{4}{n}),'fro')^2/2;
end
terms(2) = param.lambda*sum(abs(S),'all');
terms(3) = param.gamma*sum(abs(Z),'all');
temp = Y-L-S-Lam{1};
terms(4) = param.beta_1*norm(temp(mask),'fro')^2/2;
terms(5) = param.beta_2*norm(T2V(tmprod(W,D,1)-Z-Lam{2}),'fro')^2/2;
terms(6) = param.beta_3*norm(T2V(S-W-Lam{3}), 'fro')^2/2;
obj_val = sum(terms);
end