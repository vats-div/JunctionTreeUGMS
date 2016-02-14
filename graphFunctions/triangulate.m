function [cliques, G] = triangulate(G,order,ImpNodes,varargin)

% [cliques, G] = triangulate(G, order)
% 
% cliques{i} is the i'th maximal complete subgraph of the triangulated graph.
% fill_ins(i,j) = 1 iff we add a fill-in arc between i and j.
%
% To find the maximal cliques, we save each induced cluster (created by adding connecting
% neighbors) that is not a subset of any previously saved cluster. (A cluster is a complete,
% but not necessarily maximal, set of nodes.)
%
% From Kevin Murphy's BNT toolbox
% Modified by Divyanshu Vats

n = length(G);
eliminated = zeros(1,n);

if nargin > 3
    cliques = varargin{1};
    beg = length(cliques) + 1;
else
    cliques = {};
    beg = 1;
end

for i=beg:n
  u = order(i);
  U = find(~eliminated); % uneliminated
  nodes = myintersect(neighbors(G,u), U); % look up neighbors in the partially filled-in graph
  nodes = myunion(nodes, u); % the clique will always contain at least u
  G(nodes,nodes) = 1; % make them all connected to each other
  G = setdiag(G,0);  
  eliminated(u) = 1;
  
  exclude = 0;
  for c=1:length(cliques)
    if mysubset(nodes,cliques{c}) % not maximal
      exclude = 1;
      break;
    end
  end
  if ~exclude
    cnum = length(cliques)+1;
    cliques{cnum} = nodes;
  end
end

for k = 1:length(cliques)
    cliques{k} = ImpNodes(cliques{k});
end
