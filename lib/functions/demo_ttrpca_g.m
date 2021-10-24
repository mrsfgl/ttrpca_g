
dataname = 'COIL';
miss_levels = 0;%[.8:-.2:0];
noise_levels = 0;%[.5:-.1:0];
gross_levels = [0.05:0.15:0.5];
algs = {'ttrpca','ttrpca_ng'};

sz = 10*ones(1,4);
r = [4]'*ones(1,3);
lambda_list = 10.^[-1.6:.1:-.8];
alpha_list = 10.^[0];
theta_list = 10.^[-1:.1:-.5];
[X, Y, Z] = meshgrid([0.004:0.001:0.01],[0.2:.1:1],[.02:.01:.04]);
beta_list = [X(:), Y(:), Z(:)];

n_ranks = size(r,1);
n_miss = length(gross_levels);
n_exps = length(lambda_list);
results = cell(n_exps, n_ranks, n_miss);
parfor i=1:n_exps
    for i_r = 1:n_ranks
        for i_m = 1:n_miss
            results{i,i_r,i_m} = test_parameters(dataname, algs, alpha_list,...
                theta_list, lambda_list(i), beta_list,'w_ket', true,'rnd_seed', i,...
                'gross_noise', gross_levels(i_m),'miss_level', miss_levels(1),...
                'sizes', sz, 'ranks', r(i_r,:));
        end
    end
end
save(['param_',dataname,'_gross_16_lin_theta_fin.mat'],'results')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Older Experiments %%
% % % err = ones(n_miss,n_ranks,4);
% % % psnrval = ones(n_miss,n_ranks,4);
% % % ssimval = ones(n_miss,n_ranks,4);
% % % for j = 1:n_miss
% % %     for k=1:n_ranks
% % %         e_horpca = ones(n_exps, size(beta_list,1));
% % %         ps_horpca = ones(n_exps, size(beta_list,1));
% % %         ssi_horpca = ones(n_exps, size(beta_list,1));
% % %         for i=1:n_exps
% % %             e_horpca(i,:) = results{i,k,j}.err(1,1,1,1,:);
% % %             ps_horpca(i,:) = results{i,k,j}.snr(1,1,1,1,:);
% % %             ssi_horpca(i,:) = results{i,k,j}.ssim(1,1,1,1,:);
% % %         end
% % %         [err(j,k,1),in] = min(e_horpca,[],'all','linear');
% % %         psnrval(j,k,1) = max(ps_horpca,[],'all','linear');
% % %         ssimval(j,k,1) = max(ssi_horpca,[],'all','linear');
% % % %         [l_i, be_i] = ind2sub(size(e_horpca), in);
% % % %         lam_h = lambda_list(l_i);
% % % %         bet_h = beta_list(be_i,:);
% % % 
% % %         e_ttrpca = ones(n_exps, length(alpha_list), size(beta_list,1));
% % %         ps_ttrpca = ones(n_exps, length(alpha_list), size(beta_list,1));
% % %         ssi_ttrpca = ones(n_exps, length(alpha_list), size(beta_list,1));
% % %         for i=1:n_exps
% % %             e_ttrpca(i,:,:) = results{i,k,j}.err(2,:,1,1,:);
% % %             ps_ttrpca(i,:,:) = results{i,k,j}.snr(2,:,1,1,:);
% % %             ssi_ttrpca(i,:,:) = results{i,k,j}.ssim(2,:,1,1,:);
% % %         end
% % %         [err(j,k,2),in] = min(e_ttrpca,[],'all','linear');
% % %         psnrval(j,k,2) = max(ps_ttrpca,[],'all','linear');
% % %         ssimval(j,k,2) = max(ssi_ttrpca,[],'all','linear');
% % % %         [l_i, al_i, be_i] = ind2sub(size(e_ttrpca), in);
% % % %         lam_t = lambda_list(l_i);
% % % %         al_t = alpha_list(al_i);
% % % %         bet_t = beta_list(be_i,:);
% % % 
% % %         e_ttrpca_g = ones(n_exps, length(alpha_list), length(theta_list),  size(beta_list,1));
% % %         ps_ttrpca_g = ones(n_exps, length(alpha_list), length(theta_list),  size(beta_list,1));
% % %         ssi_ttrpca_g = ones(n_exps, length(alpha_list), length(theta_list),  size(beta_list,1));
% % %         for i=1:n_exps
% % %             e_ttrpca_g(i,:,:,:) = results{i,k,j}.err(3,:,:,1,:);
% % %             ps_ttrpca_g(i,:,:,:) = results{i,k,j}.snr(3,:,:,1,:);
% % %             ssi_ttrpca_g(i,:,:,:) = results{i,k,j}.ssim(3,:,:,1,:);
% % %         end
% % %         [err(j,k,3),in] = min(e_ttrpca_g,[],'all','linear');
% % %         psnrval(j,k,3) = max(ps_ttrpca_g,[],'all','linear');
% % %         ssimval(j,k,3) = max(ssi_ttrpca_g,[],'all','linear');
% % % %         [l_i, al_i, th_i, be_i] = ind2sub(size(e_ttrpca_g), in);
% % % %         lam_tg = lambda_list(l_i);
% % % %         al_tg = alpha_list(al_i);
% % % %         the_tg = theta_list(th_i);
% % % %         bet_tg = beta_list(be_i,:);
% % %         e_ttrpca_ng = ones(n_exps, length(alpha_list), length(theta_list), size(beta_list,1));
% % %         ps_ttrpca_ng = ones(n_exps, length(alpha_list), length(theta_list), size(beta_list,1));
% % %         ssi_ttrpca_ng = ones(n_exps, length(alpha_list), length(theta_list), size(beta_list,1));
% % %         for i=1:n_exps
% % %             e_ttrpca_ng(i,:,:,:) = results{i,k,j}.err(4,:,:,1,:);
% % %             ps_ttrpca_ng(i,:,:,:) = results{i,k,j}.snr(4,:,:,1,:);
% % %             ssi_ttrpca_ng(i,:,:,:) = results{i,k,j}.ssim(4,:,:,1,:);
% % %         end
% % %         [err(j,k,4),in] = min(e_ttrpca_ng,[],'all','linear');
% % %         psnrval(j,k,4) = max(ps_ttrpca_ng,[],'all','linear');
% % %         ssimval(j,k,4) = max(ssi_ttrpca_ng,[],'all','linear');
% % % %         [l_i, al_i, th_i, be_i] = ind2sub(size(e_ttrpca_ng), in);
% % % %         lam_tng = lambda_list(l_i);
% % % %         al_tng = alpha_list(al_i);
% % % %         the_tng = theta_list(th_i);
% % % %         bet_tng = beta_list(be_i,:);
% % %     end
% % % end
% % % 
% % % for i = 1:4
% % %     figure
% % %     im = imagesc(r([1,9],1), miss_levels([1,7]),1-err(:,:,i), [0,1]);
% % %     title(algs{i})
% % %     ylabel('Missing Level %%')
% % %     xlabel('Rank')
% % %     im.Interpolation = 'bilinear';
% % %     colormap(gray)
% % % end
% miss_levels = [.95:-.15:.05];
% algs = {'horpca', 'ttrpca', 'ttrpca_g', 'ttrpca_ng'};
% r = (1:9)'*ones(1,3);
% 
% for i =1:length(sz)-1
%     delta(i) = min(prod(sz(1:i)),prod(sz(i+1:end)));
% end
% delta = delta/sum(delta);
% for i=1:length(miss_levels)
%     for j=1:size(r,1)
%         data = load_data(dataname, 'miss_level', miss_levels(i),...
%         'sizes',sz,'ranks', r(j,:));
%     
%         L = horpca(data.Y, 'ind_miss', data.ind_miss,...
%                 'beta', bet_h, 'lambda', lam_h);
%         err(i,j,1) = get_SNR(L, data.Y0);
%         alpha = al_t.*delta;
%         L = ttrpca(data.Y, 'ind_miss', data.ind_miss,...
%                 'beta', [bet_t(1),bet_t(2)*alpha], 'lambda',...
%                 lam_t, 'alpha', alpha);
%         err(i,j,2) = get_SNR(L, data.Y0);
%         alpha = al_tg.*delta;
%         theta = the_tg.*delta;
%         L = ttrpca_g(data.Y, 'ind_miss', data.ind_miss,...
%                 'beta', [bet_tg(1),bet_tg(2)*alpha,bet_tg(2)*alpha], 'lambda',...
%                 lam_tg, 'alpha', alpha, 'theta', theta);
%         err(i,j,3) = get_SNR(L, data.Y0);
%         alpha = al_tng.*delta;
%         theta = the_tng.*ones(1,length(sz)-1);
%         L = ttrpca_ng(data.Y, 'ind_miss', data.ind_miss,...
%                 'beta',[bet_tng(1),bet_tng(2)*alpha,bet_tng(2)*alpha], 'lambda',...
%                 lam_tng, 'alpha', alpha, 'theta', theta);
%         err(i,j,4) = get_SNR(L, data.Y0);
%     end
% end

% save('param_test_syn.mat','results')


% dataname = 'synthetic_tt';
% miss_levels = 0;%[.8:-.2:0];
% noise_levels = 0;%[.5:-.1:0];
% gross_levels = 0.1;
% algs = {'horpca', 'ttrpca', 'ttrpca_g', 'ttrpca_ng'};
% 
% results = test_noise_miss(dataname, algs, noise_list, miss_list, gross_list);
% 
% sz = 10*ones(1,4);
% r = 4*ones(1,3);
% lambda_list = 10.^[-.8:.2:0];%10.^[-2,-1];%[.8:-.2:0];
% alpha_list = 10.^[0];%[.5:-.1:0];
% theta_list = 10.^[-7:.5:-5];
% [X, Y] = meshgrid([0.03],[1:.1:1.3]);
% beta_list = [X(:), Y(:)];
% 
% n_ranks = size(r,1);
% n_miss = length(miss_levels);
% n_exps = length(lambda_list);
% results = cell(n_exps, n_ranks, n_miss);
% parfor i=1:n_exps
%     for i_r = 1:n_ranks
%         for i_m = 1:n_miss
%             results{i,i_r,i_m} = test_parameters(dataname, algs, alpha_list,...
%                 theta_list, lambda_list(i), beta_list,'w_ket', false,'rnd_seed', i,...
%                 'gross_noise', gross_levels,'miss_level', miss_levels(i_m),...
%                 'sizes', sz, 'ranks', r(i_r,:));
%         end
%     end
% end
% 
% save('param_test_syn_2.mat','results')


% dataname = 'COIL';
% miss_levels = 0;%[.8:-.2:0];
% noise_levels = 0;%[.5:-.1:0];
% gross_levels = [0.2, 0.3, 0.4];
% algs = {'ttrpca','ttrpca_ng'};
% 
% sz = 10*ones(1,4);
% r = [4]'*ones(1,3);
% lambda_list = 10.^[-2:.2:-.4];
% alpha_list = 10.^[0];
% theta_list = 10.^[-1.5:.2:-.1];
% [X, Y, Z] = meshgrid([0.004:0.001:0.01],[0.2:.1:1],[.02:.01:.06]);
% beta_list = [X(:), Y(:), Z(:)];
% 
% n_ranks = size(r,1);
% n_miss = length(gross_levels);
% n_exps = length(lambda_list);
% results = cell(n_exps, n_ranks, n_miss);
% parfor i=1:n_exps
%     for i_r = 1:n_ranks
%         for i_m = 1:n_miss
%             results{i,i_r,i_m} = test_parameters(dataname, algs, alpha_list,...
%                 theta_list, lambda_list(i), beta_list,'w_ket', true,'rnd_seed', i,...
%                 'gross_noise', gross_levels(i_m),'miss_level', miss_levels(1),...
%                 'sizes', sz, 'ranks', r(i_r,:));
%         end
%     end
% end
% save(['param_',dataname,'_gross_16_lin_theta.mat'],'results')