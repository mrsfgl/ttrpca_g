function [A] = Lunfold(A, last_mode)
% [A] = Lunfold(A, varargin)
% Left unfolding operator.
if nargin ==2
    A = reshape(A,[],size(A,last_mode));
else
    A = reshape(A,[],size(A,ndims(A)));
end
end

