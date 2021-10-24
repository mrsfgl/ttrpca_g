function [Lx, S, obj_val, times] = ttrpca(Y, varargin)
% TTRPCA Extrats low rank and sparse tensors using 
% tensor train robust principal component analysis.
%
% [L, S, obj_val] = TTRPCA(Y)
% 


sz = size(Y);
% Parse inputs
param = inputParser;
param.addParameter('max_iter', 100);
param.addParameter('err_tol', 1e-2);
param.addParameter('alpha', [1,1,1]);
param.addParameter('lambda', 1/sqrt(max(sz)));
param.addParameter('beta', 1/(5*std(Y(:))).*ones(1,length(sz)));
param.addParameter('ind_miss', []);
param.parse(varargin{:});

max_iter = param.Results.max_iter;
err_tol = param.Results.err_tol;
alpha = param.Results.alpha;
lambda = param.Results.lambda;
beta_1 = param.Results.beta(1);
beta_2 = param.Results.beta(2:length(sz));
ind_miss = param.Results.ind_miss;

mask_Y = ones(sz)>0;
mask_Y(ind_miss) = ~mask_Y(ind_miss);
N = ndims(Y);

L = zeros(sz);
Lx = cell(1, N-1);
for n=1:N-1
    Lx{n} = zeros(sz);
end
S = zeros(sz);
Lambda{1} = zeros(sz);
Lambda{2} = cell(1, N-1);
for n=1:N-1
    Lambda{2}{n} = zeros(sz);
end

times = [];
iter = 1;
nuc_norm = num2cell(zeros(1,N));
obj_val = compute_obj(Y,L,Lx,S,Lambda,param,nuc_norm);
while true
    %% L Update
    tstart = tic;
    T2 = zeros(size(Y));
    for n=1:N-1
        T2 = T2 + beta_2(n)*(Lx{n}+Lambda{2}{n});
    end
    T1 = Y-S-Lambda{1}; 
    L(mask_Y) = (beta_1*T1(mask_Y)+T2(mask_Y))/(beta_1+sum(beta_2));
    L(~mask_Y) = T2(~mask_Y)/sum(beta_2);
    times(iter,1) = toc(tstart);
% [lobj, terms] = compute_obj(Y,L,Lx,S,Lambda,param,nuc_norm)

    %% Lx Update
    tstart = tic;
    [Lx, nuc_norm] = soft_tt_hosvd(L, Lambda{2}, alpha, beta_2);
    times(iter,2) = toc(tstart);
% [lxobj, terms] = compute_obj(Y,L,Lx,S,Lambda,param,nuc_norm)

    %% S Update
    tstart = tic;
    Sold = S;
    S = soft_threshold(Y-L-Lambda{1}, lambda/beta_1);
    times(iter,4) = toc(tstart);
% [sobj, terms] = compute_obj(Y,L,Lx,S,Lambda,param,nuc_norm)
    
    %% Dual Updates
    tstart = tic;
    temp = Y-L-S;
    Lambda{1} = Lambda{1} - temp;
    temp = norm(temp(mask_Y))/norm(Y(mask_Y));
    for n=1:N-1
        Lambda{2}{n} = Lambda{2}{n}+(Lx{n}-L);
    end
    times(iter,5) = toc(tstart);
    
    %% Error and objective calculations
    obj_val(iter+1) = compute_obj(Y,L,Lx,S,Lambda,param,nuc_norm);
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

function [val, term] = compute_obj(Y,L,Lx,S,Lambda,param,nuc_norm)
% COMPUTE_OBJ Computes the current objective value.
%
% [val, term] = compute_obj(Y,L,Lx,S,Lambda,param,nuc_norm)

alpha = param.Results.alpha;
lambda = param.Results.lambda;
beta_1 = param.Results.beta(1);
beta_2 = param.Results.beta(2:end);
N = length(Lx);
term = zeros(1,4);
for i=1:N
    term(1) = term(1) + alpha(i)*nuc_norm{i};
    term(4) = term(4) + beta_2(i)/2*sum((L-Lx{i}-Lambda{2}{i}).^2,'all');
end
term(2) = lambda*sum(abs(S),'all');
term(3) = beta_1/2*sum((Y-S-L-Lambda{1}).^2,'all');
val = sum(term);
end