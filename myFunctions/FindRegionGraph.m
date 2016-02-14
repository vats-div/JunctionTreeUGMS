function regionGraph = FindRegionGraph(JT,l,RG)

% Find region graph to process regions up to label l
% regionGraph.regions
% regionGraph.labels ->labels{k} has indices for all regions with label k
% regionGraph.edgeMat -> directed graph over regions

if nargin == 2
    regions = JT.clusters;
    labels{1} = 1:length(regions);

    % Precompute all separators in a cell
    edgesJT = triu(JT.edges); % all edges
    [indi,indj] = find(edgesJT > 0); % find all edges in JT
    numEdges = length(indi);
    edgeMat = [];

    numReg = length(regions);
    for k = 1:numEdges
        sepTemp = myintersect(JT.clusters{indi(k)},JT.clusters{indj(k)});
        if length(sepTemp) > 1
            numReg = numReg + 1;
            regions{numReg} = sepTemp;
            edgeMat(indi(k),numReg) = 1;
            edgeMat(indj(k),numReg) = 1;
        end
    end

    labels{2} = length(labels{1})+1:length(regions);

    % always find l+1 regions
    for k = 3:l+1
        % find all possible intersections of all regions in label{k-1}
        labelsTemp = labels{k-1};
        indMoreRegions = 0;
        for a1 = 1:length(labelsTemp)-1
            for a2 = a1+1:length(labelsTemp)
                sepTemp = myintersect(regions{labelsTemp(a1)},regions{labelsTemp(a2)});
                if length(sepTemp) > 1
                    % check if sepTemp is already in the regions
                    [inReg, regionNum] = IsInRegion(sepTemp,regions,labelsTemp);
                    if inReg == 0
                        indMoreRegions = 1;
                        numReg = numReg + 1;
                        regions{numReg} = sepTemp;
                        edgeMat(labelsTemp(a1),numReg) = 1;
                        edgeMat(labelsTemp(a2),numReg) = 1;
                    else
                        edgeMat(labelsTemp(a1),regionNum) = 1;
                        edgeMat(labelsTemp(a2),regionNum) = 1;
                    end
                end
            end
        end
        if indMoreRegions == 0
            break;
        end
        labels{k} = labelsTemp(end)+1:length(regions);
    end
else
    labels = RG.labels;
    regions = RG.regions;
    edgeMat = RG.edgeMat;
    labelsTemp = labels{end};
    numReg = length(regions);
    indMoreRegions = 0;
    for a1 = 1:length(labelsTemp)-1
        for a2 = a1+1:length(labelsTemp)
            sepTemp = myintersect(regions{labelsTemp(a1)},regions{labelsTemp(a2)});
            if length(sepTemp) > 1
                % check if sepTemp is already in the regions
                [inReg, regionNum] = IsInRegion(sepTemp,regions,labelsTemp);
                if inReg == 0
                    indMoreRegions = 1;
                    numReg = numReg + 1;
                    regions{numReg} = sepTemp;
                    edgeMat(labelsTemp(a1),numReg) = 1;
                    edgeMat(labelsTemp(a2),numReg) = 1;
                else
                    edgeMat(labelsTemp(a1),regionNum) = 1;
                    edgeMat(labelsTemp(a2),regionNum) = 1;
                end
            end
        end
    end
    
    if indMoreRegions == 1
        labels{length(labels)+1} = labelsTemp(end)+1:length(regions);
    end
    
end

qmax = length(regions);
edgeMat(qmax,qmax) = 0;
edgeMat = setdiag(edgeMat,0);
regionGraph.regions = regions;
regionGraph.labels = labels;
regionGraph.edgeMat = edgeMat;

end

function [inReg, regionNum] = IsInRegion(sepTemp,regions,labelsTemp)

% if no new regions
inReg = 0;
regionNum = [];
labelsEnd = labelsTemp(end)+1;
if labelsEnd > length(regions)
    return;
end

% if new regions
for k = labelsEnd:length(regions)
    if length(sepTemp) == length(regions{k})
        if isempty(mysetdiff(sepTemp,regions{k}))
            inReg = 1;
            regionNum = k;
            return;
        end
    end
end

end

% function newreg = FindIntersections(regions,labelsTemp)

% find all possible intersections of elements in regions



%end