function [X, nuc_norm] = soft_tt_hosvd(Y, Lambda, psi, tau)
% [X] = soft_tt_hosvd(Y, Lambda, tau)
% Function that returns tensor whose singular values are soft thresholded
% using parameter tau at each mode.
N = ndims(Y);
sz = size(Y);
X = cell(1,N-1);
nuc_norm = X;
for i = 1:N-1
    [X{i}, nuc_norm{i}] = soft_moden(Y-Lambda{i}, tau(i)\psi(i), 1:i);
    X{i} = m2t(X{i}, sz, 1:i);
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