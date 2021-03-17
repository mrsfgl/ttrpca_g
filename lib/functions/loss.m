function [L, S, Nt, times] = loss(Y, param)
% [S] = low_temp_sp_dec(Y, param)
% Low rank+Temporally Smooth Sparse Decomposition

N = ndims(Y);
sz = size(Y);
mask_Y = ones(sz)>0;
mask_Y(param.ind_m) = ~mask_Y(param.ind_m);
max_iter = param.max_iter;
err_tol = param.err_tol;
alpha = param.alpha;
lambda = param.lambda;
gamma = param.gamma;
psi = param.psi;
beta_1 = param.beta_1;
beta_2 = param.beta_2;
beta_3 = param.beta_3;
L = cell(1, N);
for i=1:N
    L{i} = zeros(sz);
end
S = zeros(sz);
Nt = zeros(sz);
W = zeros(sz);
Z = zeros(sz);
D = convmtx([1,-1], size(Y,1));
D(:,end) = [];
D(end,1) = -1;
Lam1 = cell(1, N);
for i=1:N
    Lam1{i} = zeros(sz);
end
Lam2 = zeros(sz);
Lam3 = zeros(sz);

timeL = [];
timeS = []; 
timeN = [];
timeZ = [];
timeW = [];
timeDual = [];
iter = 1;
% obj_val = compute_obj(Y,L,S,N,W,Z,Lam1,Lam2,Lam3,D,lambda,gamma,alpha,beta_1,beta_2,beta_3);
while true
    %% L Update
    tstart = tic;
    temp1 = zeros(sz);
    temp2 = Y-S-Nt;
    temp1(mask_Y) = temp2(mask_Y);
    [L, ~] = soft_hosvd(temp1, Lam1, psi, 1/beta_1);
    timeL(end+1)=toc(tstart);
    %% S Update
    tstart = tic;
    temp1 = zeros(sz);
    for i=1:N
        temp1 = temp1+beta_1*(Y-L{i}-Nt+Lam1{i});
    end
    temp1(~mask_Y) = 0;
    temp2 = beta_3*(W+Lam3);
    Sold = S;
    S = soft_threshold((temp1+temp2), lambda)./(N*beta_1+beta_3);
    timeS(end+1)=toc(tstart);
    %% N update
    tstart = tic;
    Nt = zeros(size(Y));
    for i=1:N
        Nt = (beta_1/(N*beta_1+alpha)).*(Y+Lam1{i}-L{i}-S);
    end
    timeN(end+1)=toc(tstart);
    %% W Update
    tstart = tic;
    Dtemp = D'*D;
    Dtemp2 = D';
    W = (beta_3*eye(sz(1))+beta_2*Dtemp)^-1 ...
        *(beta_3*Runfold(S-Lam3)+beta_2*Dtemp2*Runfold(Z+Lam2));
    W = reshape(W, sz);
    timeW(end+1)=toc(tstart);
    
    %% Z Update
    tstart = tic;
    Z = soft_threshold(mergeTensors(D, W, 2, 1)-Lam2, gamma/beta_2);
    timeZ(end+1)=toc(tstart);
    %% Dual Updates
    tstart = tic;
    temp = 0;
    for i=1:N
        Lam1_up = Y-L{i}-S-Nt;
        temp = temp + norm(Lam1_up(mask_Y))^2;
        Lam1{i}(mask_Y) = Lam1{i}(mask_Y)+Lam1_up(mask_Y);
    end
    temp = sqrt(temp)/(sqrt(N)*norm(Y(:)));
    Lam2 = Lam2-mergeTensors(D, W, 2, 1)+Z;
    Lam3 = Lam3-S+W;
    timeDual(end+1)=toc(tstart);
    
    %% Error and objective calculations
%     obj_val(iter+1) = compute_obj(Y,L,S,Nt,W,Z,Lam1,Lam2,Lam3,D,lambda,gamma,alpha,beta_1,beta_2,beta_3);
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

%  function [val, term] = compute_obj(Y,L,S,Nt,W,Z,Lam1,Lam2,Lam3,D,lambda,gamma,alpha,beta_1,beta_2,beta_3)
% N = length(L);
% term = 0;
% lag1 = 0;
% for i=1:N
%     temp = beta_1/2*sum((Y-S-L{i}-Nt+Lam1{i}).^2,'all');
%     lag1 = lag1 + temp;
%     term = term + comp_nuclear(L{i}, i);
% end
% lag2 = (beta_2/2)*sum((mergeTensors(D, W, 2, 1)-Lam2-Z).^2,'all');
% lag3 = (beta_3/2)*sum((S-W-Lam3).^2,'all');
% term = term + lag1;
% term(2) = lambda*sum(abs(S),'all')+lag1+lag3;
% term(3) = (alpha/2)*norm(Nt(:))^2+lag1;
% term(4) = lag2+lag3;
% term(5) = gamma*sum(abs(Z),'all')+lag2;
% val = sum(term)-2*lag1-lag3-lag2;
% end