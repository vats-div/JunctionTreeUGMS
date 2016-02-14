function G = UGMS_RG(H,UGMSAlg,X,options)

% G = UGMS_RG(H,UGMSAlg,X,options)
%
% Inputs:
%   H: A graph that contains the true edges (pxp adjacency matrix)
%   UGMSAlg: An undirected graphical model selection algorithm
%   X: nxp observation matrix
%   options: various options for the UGMS Algorithm
%
% UGMSAlg(X,H,Hest,D)
% Delete edges from H and add them to G
% edges
% only estimate edges in Hest
%
% author: Divyanshu Vats

if isfield(options,'Hest')
    Hest = options.Hest;
else
    Hest = H;
end
p = size(X,2);
n = size(X,1);
jtType = options.jtType;
G = zeros(p,p);

% the maximum size of the separators
SepSize = round(n/4);

if isfield(options,'SepSize')
    SepSize = options.SepSize;
end

NodesOfInterest = 1:p; % Nodes over which edges need to be estimated
ProblemSize = zeros(p,p);

% Implements a region graph based approach

indLoop = 1;

while(indLoop == 1)
    indlabel = 1;
    display('computing junction tree');
    [JT.edges,JT.clusters] = FindJunctionTree(H+G,jtType,NodesOfInterest,SepSize);
    RG = FindRegionGraph(JT,indlabel); % region graph of JT
    % RG.regions
    % RG.labels
    % RG.edgeMat

    FoundEdges = 0;
    while (FoundEdges == 0)
        % Estimate all edges possible over indlabel
        clear rr ch Htemp;
        Htemp = Hest;
        labelsTemp = RG.labels{indlabel};
        for k = 1:length(labelsTemp)
            rr{k} = RG.regions{labelsTemp(k)}; % store regions
            chr = children(RG.edgeMat,labelsTemp(k)); % all children of the region
            for ch = chr
                Htemp(RG.regions{ch},RG.regions{ch}) = 0; % do not estimate edges over children
            end
        end
        % Process each region that contains an edge to compute
        sf = @(clus) nnz(Htemp(clus,clus)); % checks if any edge in each H(cc,cc)
        checkEdge = cellfun(sf,rr);
        % only process regions where checkEdge = 1;
        checkEdgeInd = find(checkEdge >0); %#ok<EFIND>

        % if no edges to estimate in indlabel, go to next label
        if ~isempty(checkEdgeInd)
            FoundEdges = 1;
            for cInd = checkEdgeInd
                rtemp = rr{cInd};
                HccEst = spalloc(p,p,p);
                HccEst(rtemp,rtemp) = Htemp(rtemp,rtemp); %#ok<*SPRIX> % estimate these edges
                % find the nodes to apply the algorithm to
                % ancestors of RG.region{labelsTemp(cInd)}
                options.Nodes = FindAncestorsInRegionGraph(RG.edgeMat,labelsTemp(cInd),indlabel,RG.regions);
                ProblemSize(HccEst == 1) = length(options.Nodes);
                Lest = UGMSAlg(X,H+G,HccEst,options);
                G = G + full(Lest); % add all edges in Lest to G
                H = H - HccEst; % delete all edges HccEst from H
                Hest = Hest - HccEst; % delete all edges HccEst from H
            end
        else
            indlabel = indlabel + 1;
            RG = FindRegionGraph(JT,indlabel,RG); % region graph of JT
        end
        if (all(~Hest))
            % all edges have been estimated
            indLoop = 0;
        end
    end
    
    % Update ImpNodes (this reduces the complexity of updating the region
    % graph
    sf = @(clus) nnz(Hest(clus,clus)); % checks if any edge in each H(cc,cc);
    checkEdge = cellfun(sf,RG.regions);
    checkRegion = find(checkEdge >0); %#ok<EFIND>
    NodesOfInterest1 = [];
    for creg = checkRegion
        NodesOfInterest1 = [NodesOfInterest1 RG.regions{creg}];
    end
    NodesOfInterest1 = unique(NodesOfInterest1);

    % NodesofInterest
	ind_edges = find(triu(Hest) == 1);
    [ii jj] = ind2sub([p p], ind_edges);
    neighH = FindNeighborsUndirected(H+G,unique([ii jj]));
    NodesOfInterest2 = unique(neighH);
    
    NodesOfInterest = NodesOfInterest1;
    if (length(NodesOfInterest1) > length(NodesOfInterest2))
        NodesOfInterest = NodesOfInterest2;
    end
end

end

function cc = FindAncestorsInRegionGraph(edgeMat,i,indlabel,regions)

indcc = i;

for l = 1:indlabel-1
    indcc = parents(edgeMat,indcc);
end

indcc = unique(indcc);

cc = [];
for j = indcc
    cc = [regions{j} cc];
end

cc = unique(cc);

end
