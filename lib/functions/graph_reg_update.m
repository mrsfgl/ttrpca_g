function G = graph_reg_update(L,Lambda,Psi)
%GRAPH_REG_UPDATE Updates the auxiliary variable for graph regularization.
%
% G = GRAPH_REG_UPDATE(L,Lambda_3,inv_Psi)
N = length(Lambda);
G = cell(1,N);
for i=1:N
    G{i} = Psi{i}*t2m(L{i}-Lambda{i}, i);
    G{i} = reshape(G{i}, size(L{i}));
end

end