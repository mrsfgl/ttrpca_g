function [Xn, label_X] = contaminate_data(X, num_miss, len_miss, noise_var)
% [Xn, label_X] = contaminate_data(X, num_miss, len_miss, noise_var)
% Generates missing values in the data.
% 
sz = size(X);
if nargin<3
    if nargin==1
        num_miss = 100;
    end
    len_miss = sz(1);
    noise_var = 0;
end
days = randperm(sz(2)*sz(3), num_miss);
weeks = floor(days/7)+1;
days = mod(days, 7)+1;
Xn = X+noise_var*randn(sz);
label_X = zeros(sz);
for i=1:num_miss
    ind_start = randi(sz(1)-len_miss+1);
    ind_sensor = randi(sz(4));
    Xn(ind_start:ind_start+len_miss, days(i), weeks(i), ind_sensor) = 0;
    label_X(ind_start:ind_start+len_miss, days(i), weeks(i), ind_sensor) = 1;
end
end