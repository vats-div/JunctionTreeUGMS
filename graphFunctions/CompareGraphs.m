
function dist = CompareGraphs(Gtrue,G,varargin)

% function dist = CompareGraphs(Gtrue,G,varargin)
% Compares graphs Gtrue and G
% varargin = [] => edit distance
% varargin = 1 => false positive rate
% varargin = 2 => true positive rate

% author: Divyanshu Vats
% Last modified: Dec. 8, 2011

% Transform G
%tt = cell(30,1);
%if size(Gtrue,1) > size(G,1)
%	Gtrue = GeneGrouping(G,tt);
%else
%	G = GeneGrouping(G,tt);
%end

Gtrue = setdiag(Gtrue,0);
G = setdiag(G,0);

NumEdges = sum(Gtrue(:))/2;

% dist = (NumEdges - #common edges) + false edges

dist1 = G - Gtrue .* G;
dist1 = full(sum(dist1(:))/2); % type I or false positive, not there but you say it is,  number of wrong edges detected

dist2 = full(Gtrue .* G);
dist2 = NumEdges - sum(dist2(:))/2; % type II or false negative, number of edges missed

if isempty(varargin)
    dist = dist1 + dist2;
    return;
end

if varargin{1} == 1
    dist = dist1 / (sum(G(:))/2) ; % false positive rate
    %dist = dist1;
end

if varargin{1} == 2
    dist =sum(Gtrue(:).*G(:))/2;
    dist = dist/NumEdges;
end

if varargin{1} > 2
    dist = dist1 + dist2;
end

end
