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

n = ndims(X);
param = inputParser;
param.addParameter('kernel', 'Euclidean')
param.addParameter('num_class', 1)
param.addParameter('mode', n)
param.parse(varargin{:})

m = param.Results.mode;
sz = size(X);
S = prod(sz(1:m));

C = param.Results.num_class;

X = reshape(X,S,[])';
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
%     dist_mat = norm(X(:))*ones(S);
%     data = reshape(X, [], S);
%     load neighbors.mat
%     for i=1:S
%         curr_zone = data(:,i);
%         if varargin{1}
%             [~,neighbor_ind ] = intersect(regions, neighbors{i});
%         else
%             neighbor_ind = 1:S;
%         end
%         dist_mat(i, neighbor_ind) = sum((curr_zone-data(:,neighbor_ind)).^2,1)./(norm(curr_zone)*sqrt(sum(data(:,neighbor_ind).^2,1)));
%     end
%     W = exp(-(triu(dist_mat)+triu(dist_mat,1)'));
%     W = W-eye(size(W));
%     vec = sort(W(:));
%     vec(vec==inf)=[];
%     W(W<.4*vec(end))=0;
%     W = W+eye(size(W));
% %     W = W>.4*vec(end);
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