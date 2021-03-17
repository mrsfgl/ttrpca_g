function [L, nuc_norm] = soft_hosvd_wgr(Y, Lam1, La, Lam4, psi, beta_1, beta_4)
% [L, nuc_norm] = soft_hosvd_wgr(Y, Lam1, La, Lam4, psi, beta_1, beta_4)
% Function that returns the update to graph regularized nuclear norm
% minimization on tensors using ADMM.
N = ndims(Y);
sz = size(Y);
L = cell(1,N);
nuc_norm = L;
for i = 1:N
    [L{i}, nuc_norm{i}] = soft_moden(beta_1*(Y+Lam1{i})+beta_4*(La{i}+Lam4{i}), psi(i), i );
    nuc_norm{i} = nuc_norm{i}/(beta_1+beta_4);
    L{i} = m2t(L{i}, sz, i)/(beta_1+beta_4);
end
end

function [X, nuc_norm] = soft_moden(T, tau, n)
[U, S, V] = svd(t2m(T, n), 'econ');
s = diag(S)-tau;
smask = s>0;
S = diag(s(smask));
nuc_norm = sum(s(smask));
X = U(:, smask)*S*V(:, smask)';
end