function [Lx, S, obj_val, times] = ttrpca_ng(Y, varargin)
% TTRPCA_NG Extrats low rank and sparse tensors using mode-n graph 
% regularized tensor train robust principal component analysis.
%
% [L, S, obj_val] = TTRPCA_nG(Y)
% 

sz = size(Y);
N = ndims(Y);
% Parse inputs
param = inputParser;
param.addParameter('max_iter', 100);
param.addParameter('err_tol', 1e-2);
param.addParameter('alpha', [1,1,1,1]);
param.addParameter('theta', [1,1,1,1]);
param.addParameter('lambda', 1/sqrt(max(sz)));
param.addParameter('beta', 1/(5*std(Y(:))).*ones(1,2*N-1));
param.addParameter('num_neighbors',5);
param.addParameter('ind_miss', []);
param.parse(varargin{:});

max_iter = param.Results.max_iter;
err_tol = param.Results.err_tol;
alpha = param.Results.alpha;
theta = param.Results.theta;
lambda = param.Results.lambda;
beta_1 = param.Results.beta(1);
beta_2 = param.Results.beta(2:N);
beta_3 = param.Results.beta(N+1:end);
ind_miss = param.Results.ind_miss;

N = ndims(Y);
mask_Y = ones(sz)>0;
mask_Y(ind_miss) = ~mask_Y(ind_miss);

L = zeros(sz);
Lx = cell(1, N-1);
for n=1:N-1
    Lx{n} = zeros(sz);
end
G = Lx;
S = zeros(sz);
Lambda{1} = zeros(sz);
Lambda{2} = cell(1, N-1);
for n=1:N-1
    Lambda{2}{n} = zeros(sz);
end
Lambda{3} = Lambda{2};
Lambda{3}{N} = zeros(sz);
Psi = cell(1,N);
inv_Psi = cell(1,N);
for n = 1:N
    Psi{n} = get_graphL(Y, ceil(log(sz(n))), 'mode', n);
    inv_Psi{n} = (2*theta(n)/beta_3(n)*Psi{n}+eye(size(Psi{n})))^-1;
end

times = [];
iter = 1;
nuc_norm = num2cell(zeros(1,N));
obj_val = compute_obj(Y,L,Lx,G,S,Lambda,Psi,param,nuc_norm);
while true
    %% L Update
    tstart = tic;
    T2 = zeros(size(Y));
    for n=1:N-1
        T2 = T2 + beta_2(n)*(Lx{n}+Lambda{2}{n})+beta_3(n)*(G{n}+Lambda{3}{n});
    end
    T1 = Y-S-Lambda{1}; 
    L(mask_Y) = (beta_1*T1(mask_Y)+T2(mask_Y))/(beta_1+sum(beta_2)+sum(beta_3));
    L(~mask_Y) = T2(~mask_Y)/(sum(beta_2)+sum(beta_3));
    times(iter,1) = toc(tstart);
% [lobj, terms] = compute_obj(Y,L,Lx,G,S,Lambda,Psi,param,nuc_norm)
    
    %% Lx Update
    tstart = tic;
    [Lx, nuc_norm] = soft_tt_hosvd(L, Lambda{2}, alpha, beta_2);
    times(iter,2) = toc(tstart);
% [lxobj, terms] = compute_obj(Y,L,Lx,G,S,Lambda,Psi,param,nuc_norm)
    
    %% G Update
    tstart = tic;
    for n=1:N
        G{n} = inv_Psi{n}*t2m(L-Lambda{3}{n}, n);
        G{n} = m2t(G{n}, sz, n);
    end
    times(iter,3) = toc(tstart);
% [gobj, terms] = compute_obj(Y,L,Lx,G,S,Lambda,Psi,param,nuc_norm)
    
    %% S Update
    tstart = tic;
    Sold = S;
    S = soft_threshold(Y-L-Lambda{1}, lambda/beta_1);
    times(iter,4) = toc(tstart);
% [sobj, terms] = compute_obj(Y,L,Lx,G,S,Lambda,Psi,param,nuc_norm)

    %% Dual Updates
    tstart = tic;
    temp = Y-L-S;
    Lambda{1} = Lambda{1} - temp;
    temp = norm(temp(mask_Y))/norm(Y(mask_Y));
    for n=1:N-1
        Lambda{2}{n} = Lambda{2}{n}+(Lx{n}-L);
        Lambda{3}{n} = Lambda{3}{n}+(G{n}-L);
    end
    times(iter,5) = toc(tstart);
    
    %% Error and objective calculations
    obj_val(iter+1) = compute_obj(Y,L,Lx,G,S,Lambda,Psi,param,nuc_norm);
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
for n=1:N-1
    T2 = T2+Lx{n};
end
Lx = T2/N;

end

function [val, term] = compute_obj(Y,L,Lx,G,S,Lambda,Psi,param,nuc_norm)
% COMPUTE_OBJ Computes the current objective value.
%
% [val, term] = compute_obj(Y,L,Lx,G,S,Lambda,Psi,param,nuc_norm)
N = ndims(Y);
alpha = param.Results.alpha;
theta = param.Results.theta;
lambda = param.Results.lambda;
beta_1 = param.Results.beta(1);
beta_2 = param.Results.beta(2:N);
beta_3 = param.Results.beta(N+1:end);
N = length(Lx);
term = zeros(1,6);
for i=1:N
    term(1) = term(1) + alpha(i)*nuc_norm{i};
    term(2) = term(2) + theta(i)*comp_gr_reg(G{i}, Psi{i}, i);
    term(5) = term(5) + beta_2(i)/2*sum((L-Lx{i}-Lambda{2}{i}).^2,'all');
    term(6) = term(6) + beta_3(i)/2*sum((L-G{i}-Lambda{3}{i}).^2,'all');
end
term(3) = lambda*sum(abs(S),'all');
term(4) = beta_1/2*sum((Y-S-L-Lambda{1}).^2,'all');
val = sum(term);
end