function [G, H, jttime, SearchComp] = UGMS_JT(H,UGMSAlg,X,options)

% UGMSAlg(X,H,Hest,D)
% Delete edges from H and add them to G
% edges

p = size(X,2);
jtType = options.jtType;
G = zeros(p,p);
indjtreeEdges = 0;
%jttime = 0;
SearchComp = 0;
% SepSize = p;
% SepSize = round(sqrt(p));
SepSize = round(p/4);
jttime =0;

ClusterRunOnce = 0;
if isfield(options,'ClusterRunOnce');
    ClusterRunOnce = options.ClusterRunOnce;
end


while (indjtreeEdges == 0)
    numTimesCluster = 0;
    indclus = 0;
    while (indclus == 0)
        % Update Junction Tree
        [JT.edges,JT.clusters] = FindJunctionTree(H+G,jtType,1:p,SepSize);
        numTimesCluster = numTimesCluster + 1;

        if length(JT.clusters) == 1
            G = spalloc(p,p,1); H = G; SearchComp = 0;
        end
        % Precompute all separators in a cell
        edgesJT = triu(JT.edges); % all edges
        [indi,indj] = find(edgesJT > 0); % find all edges in JT
        numEdges = length(indi);
        SepR = cell(numEdges,numEdges);
        Htemp = H;
        for k = 1:numEdges
            SepR{indi(k),indj(k)} = myintersect(JT.clusters{indi(k)},JT.clusters{indj(k)});
            SepR{indj(k),indi(k)} = SepR{indi(k),indj(k)};
            Htemp(SepR{indi(k),indj(k)},SepR{indi(k),indj(k)}) = 0;
        end

        if numTimesCluster > sqrt(p)
            break;
        end
        Hold = H;

        % Process each cluster that contains an edge to compute
        sf = @(clus) nnz(Htemp(clus,clus)); % checks if any edge in each H(cc,cc)
        %sf = @(clus) any(any(Htemp(clus,clus))); % checks if any edge in each H(cc,cc)
        checkEdge = cellfun(sf,JT.clusters);
        % only process clusters where checkEdge = 1;
        checkEdgeInd = find(checkEdge > 0);

        if isempty(checkEdgeInd)
            break;
        end

        for cInd = checkEdgeInd
            cc = JT.clusters{cInd}; % cluster of nodes
            HccEst = spalloc(p,p,p);
            HccEst(cc,cc) = H(cc,cc); %#ok<*SPRIX> % estimate these edges, which we update later
            Sep = [];
            num_sep = 0;
            nsJT = FindNeighborsUndirected(JT.edges,cInd);
            for ind_k = nsJT(1:end-1)
                num_sep = num_sep + 1;
                Sep{num_sep} = SepR{cInd,ind_k};
                HccEst(Sep{num_sep},Sep{num_sep}) = 0;
            end
            options.Nodes = cc;
            % length(cc)
            [Lest, sc] = UGMSAlg(X,H+G,HccEst,options);
            SearchComp = SearchComp + sc;
            G = G + full(Lest); % add all edges in Lest to G
            H = H - HccEst; % delete all edges HccEst from H
        end

        % Check if any change in H
        if (all(~(Hold-H)))
            break;
        end

        if (all(~H))
            % all edges have been estimated
            return;
        end
        if ClusterRunOnce == 1
            indclus = 1;
        end
    end

    % Apply UGMS on edges in the junction tree
    Hold = H;

    sf = @(ss) nnz(H(ss,ss)); % checks if any edge in each H(cc,cc)
    checkEdge = cellfun(sf,SepR);
    [indi,indj] = find(triu(checkEdge) > 0); % find all edges in JT
    numEdges = length(indi);

    for sInd = 1:numEdges
        ci = indi(sInd); cj = indj(sInd);
        se = SepR{ci,cj};
        HccEst = spalloc(p,p,p);
        HccEst(se,se) = H(se,se);
        cc1 = JT.clusters{ci};
        cc2 = JT.clusters{cj};
        cc = myunion(cc1,cc2);
        % Find all sep connected to cluster indi(sInd) and indj(sInd)
        nsJTi = setdiff(FindNeighborsUndirected(JT.edges,ci,1),cj);
        nsJTj = setdiff(FindNeighborsUndirected(JT.edges,cj,1),ci);
        num_sep = 0;
        Sep = [];
        for ind_k = nsJTi(1:end)
            num_sep = num_sep + 1;
            Sep{num_sep} = SepR{ci,ind_k};
            HccEst(Sep{num_sep},Sep{num_sep}) = 0;
        end
        for ind_k = nsJTj(1:end)
            num_sep = num_sep + 1;
            Sep{num_sep} = SepR{cj,ind_k};
            HccEst(Sep{num_sep},Sep{num_sep}) = 0;
        end
        if any(HccEst(:))
            options.Nodes = cc;
            [Lest, sc] = UGMSAlg(X,H+G,HccEst,options);
            %[Lest, sc] = UGMSAlg(X,H,HccEst,[],cc);
            SearchComp = SearchComp + sc;
            G = G + full(Lest); % add all edges in Lest to G
            H = H - HccEst; % delete all edges HccEst from H
            if (all(~H)); break; end
        end
    end

    if (all(~H))
        break;
    end

    % Check if any change in H
    if (all(~(Hold-H)))
        break;
    end
    if ClusterRunOnce == 1
        indjtreeEdges = 1;
    end
end

%if options.EstimateXi == 1
    % Process edges that are leftover

    sf = @(ss) nnz(H(ss,ss)); % checks if any edge in each H(cc,cc)
    checkEdge = cellfun(sf,SepR);
    [indi,indj] = find(triu(checkEdge) > 0); % find all edges in JT
    numEdges = length(indi);
    for k = 1:numEdges
        ci = indi(k); cj = indj(k);
        se = SepR{ci,cj};
        HccEst = spalloc(p,p,p);
        HccEst(se,se) = H(se,se);
        if any(HccEst(:))
            % find all clusters that contain at least one edge from HccEst
            cc = [];
            indc = 0;
            for ff = 1:length(JT.clusters)
                H = setdiag(H,0);
                G = setdiag(G,0);
                cls = JT.clusters{ff};
                HTemp = H(cls,cls) + G(cls,cls);
                CheckClus = (sum(sum(HccEst(cls,cls) .* HTemp)) > 0);
                if CheckClus == 1
                    indc = indc + 1;
                    cc(indc) = ff;
                end
            end
            % cc contains all edges over which we do estimation
            % now find all related separators on cc
            Sep = [];
            indsep = 0;
            Nodes = [];
            for k1 = 1:length(cc)
                Nodes = [Nodes JT.clusters{cc(k1)}];
                nsJT = setdiff(FindNeighborsUndirected(JT.edges,cc(k1),1),cc);
                for k2 = 1:length(nsJT)
                    if ~isempty(SepR{cc(k1),nsJT(k2)})
                        indsep = indsep + 1;
                        Sep{indsep} = SepR{cc(k1),nsJT(k2)};
                    end
                end
            end
            options.Nodes = unique(Nodes);
            [Lest, sc] = UGMSAlg(X,H+G,HccEst,options);
            SearchComp = SearchComp + sc;
            G = G + full(Lest); % add all edges in Lest to G
            H = H - HccEst; % delete all edges HccEst from H
            if (all(~H)); break; end
        end
    end
%end

end
