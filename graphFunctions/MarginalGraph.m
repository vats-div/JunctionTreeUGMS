function MG = MarginalGraph(G, nodes)
% Marginilize graph

p = length(G);
deg = sum(G); %degree of each node
V = 1:p;
V = V(deg > 0);
othernodes = mysetdiff(V,nodes);
for i = othernodes
    ne = FindNeighborsUndirected(G,i);
    G(ne,ne) = 1;
end

MG = G(nodes,nodes);
%MG = setdiag(MG,0);
% for i=beg:n
%   u = order(i);
%   U = find(~eliminated); % uneliminated
%   nodes = myintersect(neighbors(G,u), U); % look up neighbors in the partially filled-in graph
%   nodes = myunion(nodes, u); % the clique will always contain at least u
%   G(nodes,nodes) = 1; % make them all connected to each other
%   G = setdiag(G,0);  
%   eliminated(u) = 1;
%   
%   exclude = 0;
%   for c=1:length(cliques)
%     if mysubset(nodes,cliques{c}) % not maximal
%       exclude = 1;
%       break;
%     end
%   end
%   if ~exclude
%     cnum = length(cliques)+1;
%     cliques{cnum} = nodes;
%   end
% end
