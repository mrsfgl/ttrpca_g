function val=comp_nuclear(T, n)
[~, S, ~] = svd(t2m(T,n), 'econ');
val = trace(S);
end