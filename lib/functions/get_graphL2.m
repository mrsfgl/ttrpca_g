function [L, L_inv] = get_graphL2(X, theta)
%
%
L = zeros(size(X,4),size(X,4),size(X,2));
L_inv = L;
for j=1:size(X,2)
        data = reshape(X(:,j,:,:),[],size(X,4));
        dist_mat = zeros(size(data,2));
        for i=1:size(data,2)
            curr_zone = data(:,i);
            dist_mat(i,:) = sum((curr_zone-data).^2,1)./(norm(curr_zone)*sum(data.^2,1));
        end
        W = exp(-(triu(dist_mat)+triu(dist_mat,1)'));
        W = W-min(W(:));
        vec = sort(W(:));
        vec(vec==inf)=[];
        W(W<.4*vec(end))=0;
        D = diag(sum(W));
        L(:,:,j) = D-W;
        L_inv(:,:,j) = ((theta)*L(:,:,j)+eye(size(data,2)))^-1;
end
% W = normc(W);
end