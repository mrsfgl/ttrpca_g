function [TO] = mergeTensors(T1, T2, modes_1, modes_2)
% [TO] = mergeTensors(T1, T2, modes_1, modes_2)
% Merges a tensor with tensor factors, or another tensor on specified
% modes.
s1 = size(T1);                                                             % Get the size of first tensor
if max(modes_1)>length(s1)                                                      
    s1(length(s1)+1:max(modes_1))=1;                                            % If the mode length is higher than the number of modes, then probably last modes are of dimension 1. Thus we add these to the size vector.
end
N  = length(s1);                                                           % Number of modes_1 = max(length(modes), size(T1)).
if nargin==3                                                       % Two merging operations are possible, merge a Tensor with factors or a tensor with another tensor.
    if ~iscell(T2)
        error('Missing merging modes for second tensor.');
    else
        l = length(modes_1);
        TO = mergeTensors(T1, T2{1}, modes_1(1), 2);
        for i = 2:l
            modes_1(i) = modes_1(i)-sum(modes_1(1:i-1)<modes_1(i));
            TO = mergeTensors(TO, T2{i}, [N-i+3, modes_1(i)], [1,2]);
        end
    end
else
    md1= setdiff(1:N, modes_1);
    s2 = size(T2);
    if max(modes_2)>length(s2)
        s2(length(s2):max(modes_2))=1;
    end
    N2 = length(s2); md2 = setdiff(1:N2, modes_2);
    if s1(modes_1)~=s2(modes_2)
        error('Sizes of merged modes do not match!');
    end
    T2 = reshape(permute(T2, [modes_2, md2]), prod(s2(modes_2)), []);
    T1 = reshape(permute(T1, [md1, modes_1]), [], prod(s1(modes_1)));
    if ~isempty(md1) && ~isempty(md2)
        TO = reshape(T1*T2, [s1(md1), s2(md2)]);
    elseif ~isempty(md1) && isempty(md2)
        TO = reshape(T1*T2, s1(md1));
    elseif isempty(md1) && ~isempty(md2)
        TO = reshape(T1*T2, s2(md2));
    else
        TO = T1*T2;
    end
    if length(s1)==2 && s1(end)==1                                         % If the input is a column vector and the operation is not on rows, remove the excess size from the operation.
        if modes_1==1
            TO = permute(TO, [2:ndims(TO),1]);
        end
    end
end

end

