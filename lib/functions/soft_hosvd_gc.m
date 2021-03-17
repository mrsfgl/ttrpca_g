function [L, nuc_norm] = soft_hosvd_gc(Y, Lam, Sig, U, psi, beta)
% [L, nuc_norm] = soft_hosvd_gc(Y, Lam, Sig, U, psi, beta)
% Function that returns the update to nuclear norm minimization of graph
% core using ADMM.
N = ndims(Y);
sz = cellfun(@size, U, num2cell(2*ones(1,length(U))));
L = cell(1,N);
nuc_norm = L;
for i = 1:N
    tmp = beta*(Y-Lam{i});
    [tempL, nuc_norm{i}] = soft_mode_wtd(tmp, psi(i).*(Sig{i}), i );
    nuc_norm{i} = nuc_norm{i}/beta;
    L{i} = m2t(tempL, sz, i)/beta;
end
end

function [X, nuc_norm] = soft_mode_wtd(T, tau, n)
[U, S, V] = svd(t2m(T, n), 'econ');
s = diag(S)-tau;
smask = s>0;
S = diag(s(smask));
nuc_norm = sum(s(smask));
X = U(:, smask)*S*V(:, smask)';
end