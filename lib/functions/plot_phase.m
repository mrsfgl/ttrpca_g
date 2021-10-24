for i=1:n_ranks
    for j=1:n_miss
        for k=1:length(lambda_list)
            err_horpca(i,j,k) = results{1,i,j}.err(1,1,1,k,1);
            err_ttrpca(i,j,k) = results{1,i,j}.err(2,1,1,k,1);
            for ii = 1:length(theta_list)
                err_ttrpca_g(i,j,k,ii) = results{1,i,j}.err(3,1,ii,k,1);
                err_ttrpca_ng(i,j,k,ii) = results{1,i,j}.err(4,1,ii,k,1);
            end
        end
    end
end

for k=1:length(lambda_list)
%     figure
%     im = imagesc(miss_levels([1,n_miss]), r([1,n_ranks]), (1-err_horpca(:,:,k)), [0,1]);
%     title(['HORPCA w $\lambda$=', num2str(lambda_list(k))])
%     xlabel('Missing Level %%')
%     ylabel('Rank')
%     im.Interpolation = 'bilinear';
%     colormap(gray)
%     
%     figure
%     im = imagesc(miss_levels([1,n_miss]), r([1,n_ranks]), (1-err_ttrpca(:,:,k)), [0,1]);
%     title(['TTRPCA w $\lambda$=', num2str(lambda_list(k))])
%     xlabel('Missing Level %%')
%     ylabel('Rank')
%     im.Interpolation = 'bilinear';
%     colormap(gray)
    for ii = 1:length(theta_list)
        figure
        im = imagesc(miss_levels([1,n_miss]), r([1,n_ranks]),1-err_ttrpca_g(:,:,k), [0,1]);
        title(['TTRPCA-G w $\lambda$=', num2str(lambda_list(k)),' and $\theta$=', num2str(theta_list(ii))])
        xlabel('Missing Level %%')
        ylabel('Rank')
        im.Interpolation = 'bilinear';
        colormap(gray)
    
        figure
        im = imagesc(miss_levels([1,n_miss]), r([1,n_ranks]),1-err_ttrpca_ng(:,:,k), [0,1]);
        title(['TTRPCA-nG w $\lambda$=', num2str(lambda_list(k)),' and $\theta$=', num2str(theta_list(ii))])
        xlabel('Missing Level %%')
        ylabel('Rank')
        im.Interpolation = 'bilinear';
        colormap(gray)
    end
end