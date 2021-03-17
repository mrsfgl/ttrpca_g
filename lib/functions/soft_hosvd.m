function [X, nuc_norm] = soft_hosvd(Y, Lam, psi, tau)
% [X] = soft_hosvd(Y, Lam, tau)
% Function that returns tensor whose singular values are soft thresholded
% using parameter tau at each mode.
N = ndims(Y);
sz = size(Y);
X = cell(1,N);
nuc_norm = X;
for i = 1:N
    [X{i}, nuc_norm{i}] = soft_moden(Y-Lam{i}, tau*psi(i), i );
    X{i} = m2t(X{i}, sz, i);
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