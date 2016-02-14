function [G, lambdaChosen] = UGMS_NLasso(X,H,Hest,options)

% function [G, lambdaChosen] = UGMS_NLasso(X,H,Hest,options)
% options.METHOD ->     1 for EBIC and 2 CV
% Estimate edges in Hest using the graph H

n = size(X,1);
p = size(X,2);
lambdaChosen = 0;
NumRegul = 50;
options.ReturnIndex = 1;
display(['Number of edges ' num2str(sum(Hest(:)))]);

if isempty(H)
    H = ones(p,p);
    options.Nodes = 1:p;
end

if isempty(Hest)
    Hest = H;
    options.Nodes = 1:p;
end

if isfield(options,'NumRegul')
    NumRegul = options.NumRegul; %#ok<NASGU>
else
    options.NumRegul = NumRegul;
end

if ~isfield(options,'Nodes')
    options.Nodes = 1:p;
end

if length(options.Nodes) < 8
    G = CheckEdgesRobust(X,H,Hest,options.Nodes);
    lambdaChosen = 0;
    return;
end

% find all edges to estimate
ind_edges = find(Hest == 1);
[indi indj] = ind2sub([p p], ind_edges); % the edge
edgeNodes = unique([indi indj]); % only apply UGMS to these nodes
cc = options.Nodes;

ccOld = cc; clear cc;
% rearrane cc
cc = edgeNodes;
ccdiff = setdiff(ccOld',cc);
cc = [cc ; ccdiff];

Hcc = MarginalGraph(H,cc);
Hest = Hest(cc,cc);

G = spalloc(p,p,p);

for k = 1:length(edgeNodes);
    i = edgeNodes(k); % for y
    y = X(:,i);
    ne = FindNeighborsUndirected(Hcc,k,1); % index of neighbors in cc
    neEst = FindNeighborsUndirected(Hest,k,1); % index of entries to estimate in ne
    [c, ia, ib] = intersect(ne,neEst); %#ok<NASGU,ASGLU>
    lambdaPattern = zeros(length(ne),1);
    lambdaPattern(ia) = 1;
    options.lambdaPattern = lambdaPattern;
    XX = X(:,cc(ne)); % defining X    
    bhatNonZero = myLasso(y,XX,options); % this will return locations of non-zero coefficients in XX
    ChooseLocationsTemp = cc(ne);
    ChooseLocations = ChooseLocationsTemp(bhatNonZero);
    G(i,ChooseLocations) = 1;
end

G = triu(G,1) + tril(G,-1)';
G = (G == 2);
G = triu(G) + triu(G)';

end