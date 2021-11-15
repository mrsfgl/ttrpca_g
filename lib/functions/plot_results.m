function [err, psnrval, ssimval] = plot_results(results, levels, ranks, exp_type)

n_exps = size(results,1);
n_levels = length(levels);
n_ranks = size(ranks,1);

algs = results{1}.algs;
n_algs = length(algs);
alpha_list = results{1}.alpha_list;
n_alpha = length(alpha_list);
lambda_list = results{1}.lambda_list;
n_lambda = length(lambda_list);
theta_list = results{1}.theta_list;
n_theta = length(theta_list);
beta_list = results{1}.beta_list;
n_beta = size(beta_list,1);

err = ones(n_levels,n_ranks,n_algs);
psnrval = ones(n_levels,n_ranks,n_algs);
ssimval = ones(n_levels,n_ranks,n_algs);

e_horpca = ones(n_lambda, n_beta, n_exps, n_levels);
ps_horpca = ones(n_lambda, n_beta, n_exps, n_levels);
ssi_horpca = ones(n_lambda, n_beta, n_exps, n_levels);

e_ttrpca = ones(n_alpha, n_lambda, n_beta, n_exps, n_levels);
ps_ttrpca = ones(n_alpha, n_lambda, n_beta, n_exps, n_levels);
ssi_ttrpca = ones(n_alpha, n_lambda, n_beta, n_exps, n_levels);

e_ttrpca_g = ones(n_alpha, n_theta, n_lambda, n_beta, n_exps, n_levels);
ps_ttrpca_g = ones(n_alpha, n_theta, n_lambda, n_beta, n_exps, n_levels);
ssi_ttrpca_g = ones(n_alpha, n_theta, n_lambda, n_beta, n_exps, n_levels);

e_ttrpca_ng = ones(n_alpha, n_theta, n_lambda, n_beta, n_exps, n_levels);
ps_ttrpca_ng = ones(n_alpha, n_theta, n_lambda, n_beta, n_exps, n_levels);
ssi_ttrpca_ng = ones(n_alpha, n_theta, n_lambda, n_beta, n_exps, n_levels);
for j = 1:n_levels
    for k=1:n_ranks
        if any(matches(algs, 'HORPCA', 'IgnoreCase', true))
            alg_ind = find(matches(algs, 'HORPCA', 'IgnoreCase', true));
            for i=1:n_exps
                e_horpca(:,:,i,j) = results{i,k,j}.err(alg_ind,1,1,:,:);
                ps_horpca(:,:,i,j) = results{i,k,j}.snr(alg_ind,1,1,:,:);
                ssi_horpca(:,:,i,j) = results{i,k,j}.ssim(alg_ind,1,1,:,:);
            end
            [err(j,k,alg_ind),in] = min(e_horpca,[],'all','linear');
            psnrval(j,k,alg_ind) = max(ps_horpca,[],'all','linear');
            ssimval(j,k,alg_ind) = max(ssi_horpca,[],'all','linear');
%         [l_i, be_i] = ind2sub(size(e_horpca), in);
%         lam_h = lambda_list(l_i);
%         bet_h = beta_list(be_i,:);
        end

        if any(matches(algs, 'TTRPCA','IgnoreCase',true))
            % TTRPCA
            alg_ind = find(matches(algs, 'TTRPCA', 'IgnoreCase', true));
            for i=1:n_exps
                e_ttrpca(:,:,:,i,j) = results{i,k,j}.err(alg_ind,:,1,:,:);
                ps_ttrpca(:,:,:,i,j) = results{i,k,j}.snr(alg_ind,:,1,:,:);
                ssi_ttrpca(:,:,:,i,j) = results{i,k,j}.ssim(alg_ind,:,1,:,:);
            end
            m_e_ttrpca = mean(e_ttrpca,4);
            m_ps_ttrpca = mean(ps_ttrpca,4);
            m_ssi_ttrpca = mean(ssi_ttrpca,4);
            [err(j,k,alg_ind),in] = min(m_e_ttrpca ,[],'all','linear');
            psnrval(j,k,alg_ind) = max(m_ps_ttrpca,[],'all','linear');
            ssimval(j,k,alg_ind) = max(m_ssi_ttrpca,[],'all','linear');
    %         [l_i, al_i, be_i] = ind2sub(size(e_ttrpca), in);
    %         lam_t = lambda_list(l_i);
    %         al_t = alpha_list(al_i);
    %         bet_t = beta_list(be_i,:);w
        end

        if any(matches(algs, 'TTRPCA_G', 'IgnoreCase', true))
            % TTRPCA_G
            alg_ind = find(matches(algs, 'TTRPCA_G', 'IgnoreCase', true));
            for i=1:n_exps
                e_ttrpca_g(:,:,:,:,i,j) = results{i,k,j}.err(alg_ind,:,:,:,:);
                ps_ttrpca_g(:,:,:,:,i,j) = results{i,k,j}.snr(alg_ind,:,:,:,:);
                ssi_ttrpca_g(:,:,:,:,i,j) = results{i,k,j}.ssim(alg_ind,:,:,:,:);
            end
            m_e_ttrpca_g = mean(e_ttrpca_g,5);
            m_ps_ttrpca_g = mean(ps_ttrpca_g,5);
            m_ssi_ttrpca_g = mean(ssi_ttrpca_g,5);
            [err(j,k,alg_ind),in] = min(m_e_ttrpca_g ,[],'all','linear');
            psnrval(j,k,alg_ind) = max(m_ps_ttrpca_g,[],'all','linear');
            ssimval(j,k,alg_ind) = max(m_ssi_ttrpca_g,[],'all','linear');
%             [l_i, al_i, th_i, be_i] = ind2sub(size(e_ttrpca_g), in);
%             lam_tg = lambda_list(l_i);
%             al_tg = alpha_list(al_i);
%             the_tg = theta_list(th_i);
%             bet_tg = beta_list(be_i,:);
        end
        
        if any(matches(algs, 'TTRPCA_nG', 'IgnoreCase', true))
            % TTRPCA_nG
            alg_ind = find(matches(algs, 'TTRPCA_nG', 'IgnoreCase', true));
            for i=1:n_exps
                e_ttrpca_ng(:,:,:,:,i,j) = results{i,k,j}.err(alg_ind,:,:,:,:);
                ps_ttrpca_ng(:,:,:,:,i,j) = results{i,k,j}.snr(alg_ind,:,:,:,:);
                ssi_ttrpca_ng(:,:,:,:,i,j) = results{i,k,j}.ssim(alg_ind,:,:,:,:);
            end
            m_e_ttrpca_ng = mean(e_ttrpca_ng,5);
            m_ps_ttrpca_ng = mean(ps_ttrpca_ng,5);
            m_ssi_ttrpca_ng = mean(ssi_ttrpca_ng,5);
            [err(j,k,alg_ind),in] = min(m_e_ttrpca_ng ,[],'all','linear');
            psnrval(j,k,alg_ind) = max(m_ps_ttrpca_ng,[],'all','linear');
            ssimval(j,k,alg_ind) = max(m_ssi_ttrpca_ng,[],'all','linear');
%             [l_i, al_i, th_i, be_i] = ind2sub(size(e_ttrpca_ng), in);
%             lam_tng = lambda_list(l_i);
%             al_tng = alpha_list(al_i);
%             the_tng = theta_list(th_i);
%             bet_tng = beta_list(be_i,:);
        end
    end
end

if matches(exp_type, 'phase')
    for i = 1:4
        figure
        im = imagesc(r([1,9],1), levels([1,7]),1-err(:,:,i), [0,1]);
        title(algs{i})
        ylabel('Missing Level %%')
        xlabel('Rank')
        im.Interpolation = 'bilinear';
        colormap(gray)
    end
elseif matches(exp_type, 'gross')
%     figure        
%     subplot(1,3,1)
%     suptitle('Experiments Against Gross Noise')
%     bar(levels,squeeze(err));
%     xlabel('Noise Level %')
%     ylabel('RSE')
%     legend(algs,'location','southeast')
    
%     subplot(1,3,2)
%     bar(levels,squeeze(psnrval));
%     xlabel('Noise Level %')
%     ylabel('PSNR')
%     legend(algs,'location','southeast')
    
%     subplot(1,3,3)
%     bar(levels,squeeze(ssimval));
%     xlabel('Noise Level %%')
%     ylabel('SSIM')
%     legend(algs,'location','southeast')
    
    for j=1:n_levels
        mnlvl = min(min(ps_ttrpca_ng(:,:,:,:,:,j),[],'all'), min(ps_ttrpca(:,:,:,:,j),[],'all'));
        mxlvl = max(max(ps_ttrpca_ng(:,:,:,:,:,j),[],'all'), max(ps_ttrpca(:,:,:,:,j),[],'all'));
        figure,
        for i=1:size(ps_ttrpca_ng, 2)
            subplot(1,size(ps_ttrpca_ng,2),i)
            imagesc(squeeze(ps_ttrpca_ng(:,i,:,:,:,j)),[mnlvl,mxlvl])
        end
        figure, imagesc(squeeze(ps_ttrpca(:,:,1:length(unique(beta_list(:,1:2),'rows')),:,j)),[mnlvl,mxlvl])
        figure, imagesc(squeeze(ps_horpca(:,:,1:length(unique(beta_list(:,1:2),'rows')),:,j)),[mnlvl,mxlvl])
    end
end
end