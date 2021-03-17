function [Uout] = mergeFactors(U)
%mergeFactors Merges Tensor factors.
%   Uout = mergeFactors(U)
%   U     : Is a cell array of tensor factors.
if ~iscell(U)
    error('The factors need to be a cell array.')
end
l = length(U);

Uout = U{1};
if l~=1
    if size(U{2}, 1) ~= 1
        nD = ndims(Uout);
    else 
        nD = ndims(Uout)+1;
    end
end
for i=2:l
    Uout = merge2(Uout, U{i}, nD);
    if size(U{i}, 3) == 1
        nD = ndims(Uout)+1;
    else
        nD = ndims(Uout);
    end
end

end

function U = merge2(U1, U2, varargin)
% Function that merges two tensors along one mode.
% U = merge2(U1, U2, m1, m2).
% m1 : Mode at which U1 is merged to U2. Default is the last mode of size>1.
% m2 : Mode at which U2 is merged to U1. Default is the first mode.
if length(varargin)==1
    m2 = 1;
else
    m2 = varargin{2};
end
if isempty(varargin{1})
    m1 = ndims(U1);
else
    m1 = varargin{1};
end

s1 = size(U1);
if m1 > length(s1)
    s1(length(s1)+1:m1)=1;
end
s2 = size(U2);
if m2 > length(s2)
    s2(length(s2)+1:m2)=1;
end
if s1(m1)~=s2(m2)
    error('sizes do not match.') 
end
md1= setdiff(1:length(s1), m1);
U1 = reshape(permute(U1, [md1, m1]), [], size(U1,m1));
md2= setdiff(1:length(s2), m2);
U2 = reshape(permute(U2, [m2, md2]), size(U2,m2), []);

U  = reshape(U1*U2, [s1(md1),s2(md2)]);
end