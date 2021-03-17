function [A] = mMat(A, n)
%
%
sA = size(A);
A  = reshape(A, prod(sA(1:n)), []);
end