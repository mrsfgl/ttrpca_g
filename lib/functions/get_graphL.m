function [L, W] = get_graphL(X, K, varargin)
%GET_GRAPHL Computes the graph laplacian of the mode-n unfoldings.
%
%   [L, W] = GET_GRAPHL(X, K) provides the cell array of graph Laplacians
%   L and adjacencies W. K-nn binary graphs of the last mode are
%   constructed.
%   
%   [L, W] = GET_GRAPHL(X, K, 'kernel', 'Gaussian') Constructs K-nn 
%   graphs with a Gaussian kernel.
%   
%   [L, W] = GET_GRAPHL(X, K, 'mode', n) Constructs graph for the mode n
%   
n = ndims(X);
param = inputParser;
param.addParameter('kernel', 'Euclidean')
param.addParameter('num_class', 1)
param.addParameter('mode', n)
param.parse(varargin{:})

m = param.Results.mode;
S = size(X, m);

C = param.Results.num_class;

X = t2m(X, m)';
D = zeros(S,S);
W = D;

for s=1:S
    for sp=1:S
        D(s,sp) = norm(X(:,s)-X(:,sp))^2;
    end
end

ind_knn = knnsearch(X', X', 'K', K+1);
ind_knn = ind_knn(:,2:end);
map = zeros(S);
for s=1:S
    map(s, ind_knn(s, :)) = 1;
end
map = map | map';

if strcmp(param.Results.kernel,'Gaussian')
    W = zeros(S);
    W(map) = exp(-D(map)./(2*(norm(D))));
else
    if length(K)==1
        W = map;
    elseif  length(K)==2
        D=D+diag(sum(D));
        for c=1:C
            ind        = labels==c;
            ceilD      = repmat(max(mink(D(ind,ind), K(1))),sum(ind),1);
            W(ind,ind) = D(ind,ind) <= ceilD;
            floorD     = repmat(max(mink(D(~ind, ind), K(2))),S-sum(ind),1);
            W(~ind,ind)= D(~ind,ind)<= floorD;
        end
        W = W | W';
        W = W.*kron((eye(C)-2\ones(C))*2,ones(S/C));
    end
end
D = diag(sum(W));
L = D-W;
end