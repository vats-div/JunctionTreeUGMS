function FP = FalsePositiveGraph(Gtrue,Ghat,IndRate)

% Return number of false positive edges
% assume the graphs are undirected and the matrices are symmetric
% if IndRate is specified, returns the true positive rate
% author: divyanshu vats

if isempty(Ghat)
    FP = 0;
    return;
end

NumEdges = 1;
Gtrue = setdiag(Gtrue,0);
Ghat = setdiag(Ghat,0);

if nargin == 3
    if IndRate == 1
        NumEdges = sum(Ghat(:))/2; % number of edges in the estimated graph
    end
end

FP = Ghat - Gtrue .* Ghat;

% type I or false positive, not there but you say it is,  number of wrong edges detected
FP = full(sum(FP(:))/2);

if NumEdges ~= 0
    FP = FP / NumEdges;
else
    FP = 0;
end

end