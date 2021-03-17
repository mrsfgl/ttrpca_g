function val = comp_gr_reg(La, Phi, varargin)
% val = comp_gr_reg(La, Phi, varargin)
% Computes the value of graph regularization term in the objective function.
if nargin ==2
    La = Lunfold(La);
    val = trace(La*La'*Phi);
else
    m = varargin{1};
    La = t2m(La,m);
    val = trace(La*La'*Phi);
end
end