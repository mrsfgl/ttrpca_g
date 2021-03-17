function sz = get_tt_sizes(T)
% sz = get_tt_sizes(T)
% Returns TT decomposition factor sizes as a matrix wit 3 columns.
num_factors = length(T);
sz = ones(num_factors, 3);
for i = 1:num_factors
    sz(i,1:ndims(T{i})) = size(T{i});
end
end