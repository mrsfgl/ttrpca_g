function results = test_noise_miss(dataname, algs, noise_list, miss_list, gross_list)
% TEST_NOISE_MISS Tester function for varying noise and missing data levels
% using default parameters.
%
% results = TEST_NOISE_MISS(dataname, algs, noise_list, miss_list, gross_list)

length_miss = length(miss_list);
length_noise = length(noise_list);
length_gross = length(gross_list);
length_algs = length(algs);

err = zeros(length_miss,length_noise,length_gross,length_algs);
snr = zeros(length_miss,length_noise,length_gross,length_algs);
time = zeros(length_miss,length_noise,length_gross,length_algs);
for i_m = 1:length_miss
    for i_n = 1:length_noise
        for i_g = 1:length_gross
            data = load_data(dataname, miss_list(i_m), noise_list(i_n),...
                gross_list(i_g));
            
            if matches(algs, 'HORPCA', 'IgnoreCase', true)
                % HORPCA
                alg_ind = find(matches(algs, 'HORPCA', 'IgnoreCase', true));
                [L,~,~,times] = horpca(data.Y, 'ind_miss', data.ind_miss);
                [err(i_m,i_n,i_g,alg_ind), snr(i_m,i_n,i_g,alg_ind)] = get_SNR(L, data.Y0);
                time(i_m,i_n,i_g,alg_ind) = sum(times(:));
            end
            if matches(algs, 'TTRPCA_G', 'IgnoreCase', true)
                % TTRPCA_G
                alg_ind = find(matches(algs, 'TTRPCA_G', 'IgnoreCase', true));
                [L,~,obj_g,times] = ttrpca_g(data.Y, 'ind_miss', data.ind_miss);
                [err(i_m,i_n,i_g,alg_ind), snr(i_m,i_n,i_g,alg_ind)] = get_SNR(L, data.Y0);
                time(i_m,i_n,i_g,alg_ind) = sum(times(:));
            end
            if matches(algs, 'TTRPCA_nG', 'IgnoreCase', true)
                % TTRPCA_G
                alg_ind = find(matches(algs, 'TTRPCA_nG', 'IgnoreCase', true));
                [L,~,obj_ng,times] = ttrpca_ng(data.Y,'ind_miss', data.ind_miss);
                [err(i_m,i_n,i_g,alg_ind), snr(i_m,i_n,i_g,alg_ind)] = get_SNR(L, data.Y0);
                time(i_m,i_n,i_g,alg_ind) = sum(times(:));
            end
        end
    end
end

results.miss_list = miss_list;
results.noise_list = noise_list;
results.gross_list = gross_list;
results.algs = algs;
results.err = err;
results.snr = snr;
results.time = time;
end