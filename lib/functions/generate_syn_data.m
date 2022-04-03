function data = generate_syn_data(dim, ranks, varargin)
%% synthetic data generation utility.
%
% data = generate_syn_data(dim, ranks, varargin)
% Inputs: 
% ---
%   dim: Sizes of the output tensor.
%   ranks: Ranks of the tensor.
%   varargin: Whether the inherent factorization is TT- or Tucker-based
% Returns:
% ---
%   data: Tensor with size `dim`.

if ~isempty(varargin)
    dec_type = varargin{1};
else
    dec_type = 'TT';
end

N = length(dim);
if contains(dec_type, 'TT', 'IgnoreCase', true)
    U = cell(1,N);
    U{1} = orth(randn(1*dim(1),ranks(1)));
    U{1} = reshape(U{1},1,dim(1),ranks(1));
    for n=2:N-1
        U{n} = orth(randn(ranks(n-1)*dim(n),ranks(n)));
        U{n} = reshape(U{n},[],dim(n),ranks(n));
    end
    U{N} = randn(ranks(N-1),dim(N));
    data = mergeFactors(U);
    if size(data,1) ==1 && dim(1)~=1
        data = reshape(data,dim);
    end
elseif contains(dec_type, 'Tucker', 'IgnoreCase', true)
    U = cell(1,N);
    C = randn(ranks);
    for n=1:N
        U{n} = orth(randn(dim(n),ranks(n)));
    end
    data = tmprod(C,U,1:N);
else
    error('Decomposition type is not known or not supported.')
end

end