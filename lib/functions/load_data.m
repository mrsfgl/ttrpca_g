function [data] = load_data(dataname, varargin)
% LOAD_DATA Data loader for robust PCA.
%
% data = LOAD_DATA(dataname, noise_level, miss_level, gross_noise) 
%
% Inputs: 
% ---
%   dataname: ID of the data used in the experiments.
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


%% Parse inputs
params = inputParser;
params.addParameter('noise_level', 40, @isnumeric);
params.addParameter('miss_level', 0, @isnumeric);
params.addParameter('gross_noise', 0, @isnumeric);
params.addParameter('w_ket', false, @islogical);
params.addParameter('random_seed', randi(10^3));
params.addOptional('sizes', [10,10,10,10]);
params.addOptional('ranks', [4,4,4]);
params.parse(varargin{:});

noise_level = params.Results.noise_level;
miss_level = params.Results.miss_level;
gross_noise = params.Results.gross_noise;
w_ket = params.Results.w_ket;
random_seed = params.Results.random_seed;
sizes = params.Results.sizes;
ranks = params.Results.ranks;

rng(random_seed);
switch dataname
    case 'synthetic_tt'
        data.Y0 = generate_syn_data(sizes,ranks);
        
        data.Y = data.Y0;
        data.Y(:) = awgn(data.Y(:), noise_level, 'measured');
        
        if gross_noise>0
            mask = randperm(prod(sizes), round(prod(sizes)*gross_noise));
        else
            mask =[];
        end
        data.Y(mask) = rand([length(mask),1]);
        
        ind_miss = randperm(prod(sizes), round(prod(sizes)*miss_level));
        data.Y(ind_miss) = 0;
        data.ind_miss = ind_miss;
        
    case 'synthetic_tucker'
        data.Y0 = generate_syn_data(sizes, ranks,'tucker');
        
        
        data.Y = data.Y0;
        data.Y(:) = awgn(data.Y(:), noise_level, 'measured');
        
        if gross_noise>0
            mask = randperm(prod(sizes), round(prod(sizes)*gross_noise));
        else
            mask =[];
        end
        data.Y(mask) = rand([length(mask),1]);
        
        ind_miss = randperm(prod(sizes), round(prod(sizes)*miss_level));
        data.Y(ind_miss) = 0;
        data.ind_miss = ind_miss;
        
    case 'COIL'
        load('data/CoilData.mat', 'Data')

        Data = Data/max(Data,[],'all');
        Data = Data(:,:,1:2:72,1:40);
        sz = size(Data);
        I = zeros(16,16,sz(3),sz(4));
        for i=1:sz(3)
            for j=1:sz(4)
                I(:,:,i,j) = imresize(Data(:,:,i,j),.25);
            end
        end
        Data = I;
        if gross_noise>0
            mask = randperm(numel(I), round(numel(I)*gross_noise));
        else
            mask =[];
        end
        I(mask) = rand([length(mask),1]);
        
        ind_miss = randperm(numel(I), round(numel(I)*miss_level));
        data.Y0 = Data;
        data.Y = I;
        data.Y(ind_miss) = 0;
        data.ind_miss = ind_miss;
        
    case {'lena', 'peppers', 'baboon'} 
        I = imread(['data/images/',dataname,'.png']);
        
        I = double(imresize(I,.5))/255;
        if w_ket
            data.Y0 = ket_aug(I);
        else
            data.Y0 = I;
        end
        for i=1:3
            I(:,:,i) = awgn(I(:,:,i),noise_level,'measured');
        end
        if gross_noise>0
            mask = randperm(numel(I), round(numel(I)*gross_noise));
        else
            mask =[];
        end
        I(mask) = rand([length(mask),1]);
        ind_miss = randperm(numel(I), round(numel(I)*miss_level));
        I(ind_miss) = 0;
        if w_ket
            data.Y = ket_aug(I);
        else
            data.Y = I;
        end
        data.ind_miss = ind_miss;
        
    otherwise
        error('Unknown data.')
end


end

function I = ket_aug(I)

I = reshape(I,[4*ones(1,8),3]);
P = reshape(1:8,4,2)'; 
I = permute(I,[P(:)',9]);
I = reshape(I, [16*ones(1,4),3]);
end