% function G = construct_k_jtree(p,a,NumEdges,maxsep,maxdeg)
%
% returns a sparse inverse covariance matrix
% contact: Divyanshu Vats, dvats@ima.umn.edu

function G = construct_k_jtree(p,a,maxsep,maxclus)

% initialization
G = zeros(p,p);
if nargin == 3
 	maxclus = maxsep + 1;
end

rho = 1/maxsep-0.001;

rn = 1 + ceil((maxclus-1) * rand(1));

% create k-clique first
G(1:rn,1:rn) = double(rand(rn,rn) > 0.5) * rand(1) * rho * randomSign(1);

% Maintain a list of cliques in the graph
clique{1} = (1:rn);
num_cliq = 1;
% add nodes to the graph
edge = [];

V = 1:p;
indv = rn+1;
cliqueSize = length(clique{1});
while (indv <= p)
    % pick a clique
    %atemp = hist(cliqueInd,num_cliq);
    %atempind = find(atemp < 2);
    %c_num = atempind(floor(length(atempind) * rand(1)) + 1);
    c_num = (floor(num_cliq * rand(1)) + 1);
    cc = clique{c_num}; lencc = cliqueSize(c_num);
    spsize = geornd(a)+1;
    if spsize > min(maxsep,lencc-1)
        %spsize = min(maxsep,lencc-1);
        spsize = 0;
    end
    %spsize
    % choose nodes in cc to connect to
    cccopy = cc;
    lcccopy = length(cccopy);
    sep = zeros(1,spsize);
    for d = 1:spsize
        rndtemp = ceil(lcccopy * rand(1));
        sep(d) = cccopy(rndtemp);
        cccopy(cccopy == sep(d)) = [];
        lcccopy = length(cccopy);
    end
    
    % choose nodes that will connect to sep
    NextNumClus = (ceil((maxclus - spsize) * rand(1)));
    ccNext = V(indv : min(indv + NextNumClus,p));
    %ccNext = V(indv);
    %indv = indv + 1;
    indv = indv+NextNumClus+1;
    
    % add edges between ccNext and sep
    G(ccNext,sep) = rand(1) * rho * randomSign(1); G(sep,ccNext) = G(ccNext,sep)';
    G(ccNext,ccNext) = (rand(length(ccNext)) > 0.5) * rand(1) * rho * randomSign(1);
    
    % create another clique
    num_cliq = num_cliq + 1;
    clique{num_cliq} = [sep ccNext ccNext ccNext ccNext];
    %clique{num_cliq} = [sep ccNext];
    cliqueSize(num_cliq) = length(sep) + length(ccNext);
end

% remove edges over separators

[edgesJT,JT.clusters] = FindJunctionTree(G,2,p);
[indi,indj] = find(edgesJT > 1); % find all edges in JT
num_edges = length(indi);

edgesOverSep = [];
for k = 1:num_edges
    ss = myintersect(JT.clusters{indi(k)},JT.clusters{indj(k)});
    G(ss,ss) = 0;
    GG = spalloc(p,p,p);
    GG(ss,ss) = 1;
    GG = setdiag(GG,0);
    edgesOverSep = [edgesOverSep ; find(triu(GG) > 0)];
end

G = triu(G) + triu(G)';
G = setdiag(G,0);

[indi indj] = find(triu(G) ~= 0);

for k = 1:length(indi);
    nei = FindNeighborsUndirected(G ~= 0,indi(k),1);
    nej = FindNeighborsUndirected(G ~= 0,indj(k),1);
    if (length(nei) > maxsep || length(nej) > maxsep)
        G(indi(k),indj(k)) = 0.995/(max(length(nei),length(nej)));
        G(indj(k),indi(k)) = 0.995/(max(length(nei),length(nej)));
    end
end

rndsgn = randomSign(p,p);
rndsgn = triu(rndsgn) + triu(rndsgn)';
G = rndsgn .* G;

G = setdiag(G,1);

% edgesOverSep = unique(edgesOverSep);
% G = setdiag(G,0);
% G = triu(G) + triu(G)';
% CurrentNumEdges = sum(G(:))/2;
% 
% if CurrentNumEdges == NumEdges
%     return;
% end



% if (CurrentNumEdges < NumEdges)
%     % this means we need to add more edges
%     NumEdgesToAdd  = NumEdges - CurrentNumEdges;
%     % choose edges from edgesOverSep
%     rndIndex = randperm(length(edgesOverSep));
%     AddEdges = edgesOverSep(1:rndIndex);
%     for e = 1:min(NumEdgesToAdd,length(AddEdges))
%         [i j] = ind2sub([p p], AddEdges(e)); % the edge
%         G(i,j) = 1;
%         G(j,i) = 1;
%     end
% else
%     % delete edges
%     [indi,indj] = find(triu(G) > 0); % find all edges in JT
%     NumEdgesDelete = CurrentNumEdges - NumEdges;
%     rndIndex = randperm(CurrentNumEdges);    
%     indi = indi(rndIndex(1:NumEdgesDelete));
%     indj = indj(rndIndex(1:NumEdgesDelete));
%     for e = 1:NumEdgesDelete
%         G(indi(e),indj(e)) = 0;
%         G(indj(e),indi(e)) = 0;
%     end
% end

end
