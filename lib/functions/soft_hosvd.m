function [X, nuc_norm] = soft_hosvd(Y, Lam, psi, tau)
% [X] = soft_hosvd(Y, Lam, tau)
% Function that returns tensor whose singular values are soft thresholded
% using parameter tau at each mode.
N = length(Lam);
X = cell(1,N);
nuc_norm = X;
for i = 1:N
    [X{i}, nuc_norm{i}] = soft_moden(Y-Lam{i}, tau*psi(i), i );
    nuc_norm{i} = nuc_norm{i}*psi(i);
end
end

function [X, nuc_norm] = soft_moden(T, tau, n)

sz = size(T);
[U, S, V] = svd(t2m(T, n), 'econ');
s = diag(S)-tau;
smask = s>0;
S = diag(s(smask));
nuc_norm = sum(s(smask));
X = U(:, smask)*S*V(:, smask)';
X = m2t(X, sz, n);
end