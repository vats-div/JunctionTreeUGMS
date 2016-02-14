function ns = FindNeighborsUndirected(G, i, ind)
% NEIGHBORS Find the neighbors of i in the undirected graph G
% if ind = [], then include neighbor
% if ind is nonempty, then don't include neighbor

G(i,i) = 0;
[temp,ns] = ind2sub([length(i) length(G)],find(G(i,:) > 0));
ns = myunique(ns);

if (nargin == 2)
    ns(end+1:end+length(i)) = i;
end
