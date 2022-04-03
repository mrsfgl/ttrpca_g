function [results, data] = test_parameters(dataname, algs, alpha_list, theta_list, lambda_list, beta_list, varargin)
% TEST_PARAMETERS Tester function for sensitivity analysis on all
% parameters.
% 
% [results, data] = test_parameters(dataname, algs, alpha_list, theta_list, lambda_list, beta_list, varargin)
%
% Inputs: 
% ---
%   dataname: ID of the data used in the experiments.
%   algs: Algorithms to be compared. e.g. 'horpca', 'ttrpca', 'ttrpca_g'
%   alpha_list: ADMM parameters
%   theta_list: ADMM parameters
%   lambda_list: ADMM parameters
%   beta_list: ADMM parameters
% Optional Parameters:
% ---
%   sizes: array of sizes to initialize data. Necessary for synthetic
%       data. default: [10,10,10,10]
%   ranks: array of ranks to initialize data. Necessary for synthetic data.
%       default: [4,4,4]
%   noise_level: AWGN SNR. default: 40 dB
%   miss_level: Pctg of unobserved entries. default: 0
%   gross_noise: Pctg of gross noise. default: 0
%   w_ket: Apply ket augmentation or not. default: false
%   rnd_seed: Random seed for different initialization. default: random
%       integer
% 
% Returns:
% ---
%   results: Struct with output statistics and corresponding parameter
%       list.



params = inputParser;
params.addOptional('sizes', [10,10,10,10]);
params.addOptional('ranks', [4,4,4]);
params.addParameter('noise_level', 40, @isnumeric);
params.addParameter('miss_level', 0, @isnumeric);
params.addParameter('gross_noise', 0.15, @isnumeric);
params.addParameter('w_ket', false, @islogical);
params.addParameter('rnd_seed', randi(10^3));
params.parse(varargin{:});

sizes = params.Results.sizes;
ranks = params.Results.ranks;
noise_var = params.Results.noise_level;
miss_pct = params.Results.miss_level;
gross_pct = params.Results.gross_noise;
w_ket = params.Results.w_ket;
rnd_seed = params.Results.rnd_seed;

data = load_data(dataname, 'noise_level', noise_var, ...
    'miss_level', miss_pct, 'gross_noise', gross_pct,'w_ket', w_ket,...
    'random_seed',rnd_seed,'sizes',sizes,'ranks', ranks);
[results.raw_err, results.raw_snr, results.raw_ssim] = get_SNR(data.Y, data.Y0);
results.miss = miss_pct;
results.ranks = ranks;
results.size = sizes;
results.name = 'dataname';

sz = size(data.Y);

length_alpha = length(alpha_list);
length_theta = length(theta_list);
length_lambda = length(lambda_list);
length_beta = size(beta_list, 1);
num_algs = length(algs);

err = zeros(num_algs, length_alpha, length_theta, length_lambda, length_beta);
snr_val = zeros(num_algs, length_alpha, length_theta, length_lambda, length_beta);
ssim_val = zeros(num_algs, length_alpha, length_theta, length_lambda, length_beta);
time = zeros(num_algs, length_alpha, length_theta, length_lambda, length_beta);

for i =1:ndims(data.Y)
    delta(i) = min(prod(sz(1:i)),prod(sz(i+1:end)));
    delta_2(i) = sz(i);
end
delta = delta/sum(delta);
delta_2 = delta_2./sum(delta_2);

for i_l = 1:length_lambda
    lambda = lambda_list(i_l);
    for i_b = 1:length_beta
        beta = beta_list(i_b,:);
        if any(matches(algs, 'HORPCA', 'IgnoreCase', true))
            % HORPCA
            alg_ind = find(matches(algs, 'HORPCA', 'IgnoreCase', true));
            algs{alg_ind} = 'HoRPCA';
            [L,~,~,times] = horpca(data.Y, 'ind_miss', data.ind_miss,...
                'beta', beta, 'lambda', lambda);
            [err(alg_ind,:,:,i_l,i_b), snr_val(alg_ind,:,:,i_l,i_b), ssim_val(alg_ind,:,:,i_l,i_b)] = get_SNR(L, data.Y0);
            time(alg_ind,:,:,i_l,i_b) = sum(times(:));
        end
        for i_a = 1:length_alpha
            alpha = alpha_list(i_a).*delta(1:end-1);
            if any(matches(algs, 'TTRPCA','IgnoreCase',true))
                % TTRPCA
                alg_ind = find(matches(algs, 'TTRPCA', 'IgnoreCase', true));
                algs{alg_ind} = 'TTRPCA';
                [L,~,~,times] = ttrpca(data.Y, 'ind_miss', data.ind_miss,...
                    'beta', [beta(1),beta(2)*alpha], 'lambda', lambda, 'alpha', alpha);
                [err(alg_ind,i_a,:,i_l,i_b), snr_val(alg_ind,i_a,:,i_l,i_b), ssim_val(alg_ind,i_a,:,i_l,i_b)] = get_SNR(L, data.Y0);
                time(alg_ind,i_a,:,i_l,i_b) = sum(times(:));
            end
            for i_t = 1:length_theta
                theta = theta_list(i_t).*delta(1:end-1);
                
                if any(matches(algs, 'TTRPCA_G', 'IgnoreCase', true))
                    % TTRPCA_G
                    alg_ind = find(matches(algs, 'TTRPCA_G', 'IgnoreCase', true));
                    algs{alg_ind} = 'TTRPCA_G';
                    [L,~,~,times] = ttrpca_g(data.Y, 'ind_miss', data.ind_miss,...
                        'beta', [beta(1),beta(2)*alpha,beta(3)*theta], 'lambda', lambda, 'theta', theta, ...
                        'alpha', alpha);
                    [err(alg_ind,i_a,i_t,i_l,i_b), snr_val(alg_ind,i_a,i_t,i_l,i_b), ssim_val(alg_ind,i_a,i_t,i_l,i_b)] = get_SNR(L, data.Y0);
                    time(alg_ind,i_a,i_t,i_l,i_b) = sum(times(:));
                end
                theta = theta_list(i_t).*delta_2;
                if any(matches(algs, 'TTRPCA_nG', 'IgnoreCase', true))
                    % TTRPCA_nG
                    alg_ind = find(matches(algs, 'TTRPCA_nG', 'IgnoreCase', true));
                    algs{alg_ind} = 'TTRPCA_nG';
                    [L,~,~,times] = ttrpca_ng(data.Y,'ind_miss', data.ind_miss,...
                        'beta', [beta(1),beta(2)*alpha,beta(3)*theta+eps], 'lambda', lambda, 'theta', theta, ...
                        'alpha', alpha);
                    [err(alg_ind, i_a,i_t,i_l,i_b), snr_val(alg_ind, i_a,i_t,i_l,i_b), ssim_val(alg_ind,i_a,i_t,i_l,i_b)] = get_SNR(L, data.Y0);
                    time(alg_ind, i_a,i_t,i_l,i_b) = sum(times(:));
                end
            end
        end
    end
end

results.lambda_list = lambda_list;
results.alpha_list = alpha_list;
results.theta_list = theta_list;
results.beta_list = beta_list;
results.algs = algs;
results.err = err;
results.snr = snr_val;
results.ssim = ssim_val;
results.time = time;
end