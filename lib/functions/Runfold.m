function [A] = Runfold(A, first_mode)
% [A] = Runfold(A,first_mode)
% Right unfolding operator.
if nargin==2
    A = reshape(A,size(A, first_mode),[]);
else
    A = reshape(A,size(A,1),[]);
end
end

