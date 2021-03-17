function [Lx, S, obj_val] = ttrpca_g(Y, varargin)
% TTRPCAG Extrats low rank and sparse tensors using graph regularized 
% tensor train robust principal component analysis.
%
% [L, S, obj_val] = TTRPCA_G(Y)
% 


% Parse inputs
param = inputParser;
param.addParameter('max_iter', 1000);
param.addParameter('err_tol', 1e-5);
param.addParameter('alpha', [1,1,1,1]);
param.addParameter('theta', [1,1,1,1]);
param.addParameter('beta', [1e-1, 1e-1, 1e-1]);
param.parse(varargin{:});

max_iter = param.Results.max_iter;
err_tol = param.Results.err_tol;
alpha = param.Results.alpha;
theta = param.Results.theta;
beta_1 = param.Results.beta(1);
beta_2 = param.Results.beta(2);
beta_3 = param.Results.beta(3);


N = ndims(Y);
sz = size(Y);

L = zeros(sz);
Lx = cell(1, N);
for n=1:N
    Lx{n} = zeros(sz);
end
G = Lx;
S = zeros(sz);
Nt = zeros(sz);
for n = 1:N
    Phi{n} = get_graphL_G(Y, param.Results.graph_neighbors, 'mode', n);
    inv_Phi{n} = (theta(n)/beta_3*Phi{n}+eye(size(Phi)))^-1;
end
Lambda{1} = zeros(sz);
Lambda{2} = cell(1, N);
for n=1:N
    Lambda{2}{n} = zeros(sz);
end
Lambda{3} = Lambda{2};

times = [];
iter = 1;
obj_val = compute_obj(Y,Lx,G,S,Nt,Lambda,D,Phi,param);
while true
    %% L Update
    tstart = tic;
    T2 = zeros(size(Y));
    for n=1:4
        T2 = T2 + beta_2*(Lx{n}+Lambda{2}{n})+beta_3*(G{n}+Lambda{3}{n});
    end
    T1 = Y-S-Lambda{1}; 
    L(mask_Y) = (beta_1*T1(mask_Y)+T2(mask_Y))/(beta_1+4*(beta_2+beta_3));
    L(~mask_Y) = T2(~mask_Y)/(4*(beta_2+beta_3));
    times(iter,1) = toc(tstart);
    
    %% Lx Update
    tstart = tic;
    [Lx, nuc_norm] = soft_tt_hosvd(L, Lambda{2}, alpha, 1/beta_2);
    times(iter,2) = toc(tstart);
    
    %% G Update
    for n=1:N
        G{n} = inv_Phi{n}*t2m(L-Lambda{3}{n},1:n);
        G{n} = reshape(G{n}, sz);
    end
    times(iter,3) = toc(tstart);
    
    %% S Update
    tstart = tic;
    Sold = S;
    S = soft_threshold(Y-L-Lambda{1}, 1/beta_1);
    times(iter,4) = toc(tstart);
    %% Dual Updates
    tstart = tic;
    temp = Y-L-S;
    temp = norm(temp(:))/norm(Y(:));
    Lambda{1} = Lambda{1} - temp;
    for n=1:N
        Lambda{2}{n} = Lambda{2}{n}+(Lx{n}-L);
        Lambda{3}{n} = Lambda{3}{n}+(G{n}-L);
    end
    times(iter,5) = toc(tstart);
    
    %% Error and objective calculations
    obj_val(iter+1) = compute_obj(Y,L,Lx,G,S,Lambda,Phi,param,nuc_norm);
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
T2 = zeros(sz);
for n=1:N
    T2 = T2+Lx{n};
end
Lx = T2/N;

end

function [val, term] = compute_obj(Y,L,Lx,G,S,Lambda,Phi,param,nuc_norm)
% COMPUTE_OBJ Computes the current objective value.
%
% [val, term] = compute_obj(Y,L,Lx,G,S,Lambda,Phi,param,nuc_norm)
alpha = param.Results.alpha;
theta = param.Results.theta;
beta_1 = param.Results.beta(1);
beta_2 = param.Results.beta(2);
beta_3 = param.Results.beta(3);
N = length(L);
term = zeros(1,6);
for i=1:N
    term(1) = term(1) + alpha(i)*nuc_norm{i};
    term(2) = term(2) + theta(i)*comp_gr_reg(G{i}, Phi{i}, 1:i);
    term(5) = term(5) + beta_2/2*sum((L-Lx{i}-Lambda{2}{i}).^2,'all');
    term(6) = term(6) + beta_3/2*sum((L-G{i}-Lambda{3}{i}).^2,'all');
end
term(3) = sum(abs(S),'all');
term(4) = beta_1/2*sum((Y-S-L-Lambda{1}{i}).^2,'all');
val = sum(term);
end