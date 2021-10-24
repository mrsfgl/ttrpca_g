function [L, W] = get_graphL_G(X, K, varargin)
%GET_GRAPHL_G Computes the graph laplacian of the canonical unfoldings.
%
%   [L, W] = GET_GRAPHL_G(X, K) provides the cell array of graph Laplacians
%   L and adjacencies W. K-nn binary graphs are constructed.
%   
%   [L, W] = GET_GRAPHL_G(X, K, 'kernel', 'Gaussian') Constructs K-nn 
%   graphs with a Gaussian kernel.
%   
%   [L, W] = GET_GRAPHL_G(X, K, 'mode', m) Constructs graph for the modes
%   until m.
%
%   [L, W] = GET_GRAPHL(X, K, 'num_class', n) Constructs supervised graphs
%   based on the number of classes. The samples are assumed to be ordered
%   with respect to class and each class is assumed to have the same number
%   of samples.
%   

n = ndims(X);
param = inputParser;
param.addParameter('kernel', 'Euclidean')
param.addParameter('num_class', 1)
param.addParameter('mode', n-1)
param.parse(varargin{:})

m = param.Results.mode;
sz = size(X);
S = prod(sz(1:m));

C = param.Results.num_class;

X = reshape(X,S,[])';
% K = min(size(X,2)-1, K);
D = zeros(S,S);
W = D;

ind_knn = knnsearch(X', X', 'K', K+1);
ind_knn = ind_knn(:,2:end);
K = min(size(ind_knn,2), K);

for s=1:S
    for sp=1:K
        D(s,ind_knn(s,sp)) = norm(X(:,s)-X(:,ind_knn(s,sp)))^2;
    end
end
D = D+D';
map = zeros(S);
for s=1:S
    map(s, ind_knn(s, :)) = 1;
end
D(map & map') = D(map & map')/2;
map = map | map';


if strcmp(param.Results.kernel,'Gaussian')
    W = zeros(S);
    const = floor(log10(sum(var(D))));
    gamma = 1/(S*10^const);
    W(map) = exp(-D(map)./gamma);
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
% W = normc(W);
end