function [L, W] = get_graphL(X, K, varargin)
%GET_GRAPHL Computes the graph laplacian of the mode-n unfoldings.
%
%   [L, W] = GET_GRAPHL(X, K) provides the cell array of graph Laplacians
%   L and adjacencies W. K-nn binary graphs of the last mode are
%   constructed.
%
%   [L, W] = GET_GRAPHL(X, K, 'kernel', 'Gaussian') Constructs K-nn
%   graphs with a Gaussian kernel. {'Euclidean'}
%
%   [L, W] = GET_GRAPHL(X, K, 'mode', n) Constructs graph for mode n. {1}
%
%   [L, W] = GET_GRAPHL(X, K, 'num_class', n) Constructs supervised graphs
%   based on the number of classes. The samples are assumed to be ordered
%   with respect to class and each class is assumed to have the same number
%   of samples.
%
param = inputParser;
param.addParameter('kernel', 'Euclidean')
param.addParameter('num_class', 1)
param.addParameter('mode', 1)
param.parse(varargin{:})

m = param.Results.mode;
S = size(X, m);

C = param.Results.num_class;

X = t2m(X, m)';
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
    gamma = (S*10^const);
%     gamma = 2;
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

end