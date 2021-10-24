function [L, S, obj_val, times] = horpca(Y, varargin)
% HORPCA Extrats low rank and sparse tensors using higher order robust
% principal component analysis.
%
% [L, S, obj_val, times] = horpca(Y)


N = ndims(Y);

% Parse inputs
param = inputParser;
param.addParameter('max_iter', 100);
param.addParameter('ind_miss', []);
param.addParameter('err_tol', 1e-5);
param.addParameter('lambda', 10);%/sqrt(max(size(Y))));
param.addParameter('psi', ones(1,N));
param.addParameter('beta', [10,1]);%[1,1]./(5*std(Y(:))));
param.parse(varargin{:});

mask_Y = ones(size(Y))>0;
mask_Y(param.Results.ind_miss) = ~mask_Y(param.Results.ind_miss);
max_iter = param.Results.max_iter;
err_tol = param.Results.err_tol;
lambda = param.Results.lambda;
psi = param.Results.psi;
beta_1 = param.Results.beta(1);
beta_2 = param.Results.beta(2);

sz = size(Y);
Lx = cell(1, N);
for i=1:N
    Lx{i} = zeros(sz);
end
S = zeros(sz);
L = zeros(sz);

Lambda{1} = zeros(sz);
Lambda{2} = cell(1, N);
for i=1:N
    Lambda{2}{i} = zeros(sz);
end

tnorm = norm(Y(:));
iter = 1;
nuc_norm = num2cell(zeros(1,N));
obj_val = compute_obj(Y,L,Lx,S,Lambda,param.Results,nuc_norm);
while true
    %% L Update
    tstart = tic;
    T1 = Y-S-Lambda{1};
    T1(~mask_Y) = 0;
    T2 = zeros(sz);
    for n=1:N
        T2 = T2+beta_2*(Lx{n}+Lambda{2}{n});
    end
    L(mask_Y) = (beta_1*T1(mask_Y)+T2(mask_Y))/(beta_1+N*beta_2);
    L(~mask_Y) = T2(~mask_Y)/(N*beta_2);
    times(1,iter) = toc(tstart);
    
    %% Lx Update
    tstart = tic;
    [Lx, nuc_norm] = soft_hosvd(L, Lambda{2}, psi, 1/beta_2);
    times(2,iter) = toc(tstart);
    
    %% S Update
    tstart = tic;
    Sold = S;
    S = soft_threshold(Y-L-Lambda{1}, lambda/beta_1);
    times(3,iter) = toc(tstart);
    
    %% Dual Updates
    tstart = tic;
    Lambda{1} = Lambda{1}+L+S-Y;
    temp = 0;
    for i=1:N
        Lam2_up = L-Lx{i};
        temp = temp + norm(Lam2_up(:))^2;
        Lambda{2}{i} = Lambda{2}{i}-Lam2_up;
    end
    temp = sqrt(temp)/(sqrt(N)*tnorm);
    times(4,iter) = toc(tstart);
    
    %% Error and objective calculations
    obj_val(iter+1) = compute_obj(Y,L,Lx,S,Lambda,param.Results, nuc_norm);
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

end

function [val, term] = compute_obj(Y,L,Lx,S,Lambda,param,nuc_norm)

lambda = param.lambda;
beta_1 = param.beta(1);
beta_2 = param.beta(2);
N = length(Lx);
term = zeros(1,4);
term(1) = beta_1/2*sum((Y-L-S-Lambda{1}).^2,'all');
for i=1:N
    term(2) = term(2) + nuc_norm{i};
    term(4) = term(4) + beta_2/2*sum((L-Lx{i}-Lambda{2}{i}).^2,'all');
end
term(3) = lambda*sum(abs(S),'all');
val = sum(term);
end