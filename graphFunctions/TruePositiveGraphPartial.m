function TP = TruePositiveGraphPartial(Gtrue,Ghat,IndRate,p1)

% Return number of true positive edges,
% assumes the graphs are undirected and the matrices are symmetric
% if IndRate is specified, returns the true positive rate
% author: divyanshu vats

if isempty(Ghat)
    TP = 0;
    return;
end


Gtrue = Gtrue(1:p1,1:p1);
Ghat = Ghat(1:p1,1:p1);


NumEdges = 1;
Gtrue = setdiag(Gtrue,0);
Ghat = setdiag(Ghat,0);

%if nargin == 3
    if IndRate == 1
        NumEdges = sum(Gtrue(:))/2;
    end
%end

TP =sum(Gtrue(:).*Ghat(:))/2; % number of common edges
TP = TP/NumEdges;

end