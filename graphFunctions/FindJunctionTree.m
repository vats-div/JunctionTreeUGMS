function [jtree, cliques, MG, Gtriag] = FindJunctionTree(MG,varargin)

% This is a modified version of the code in Kevin Murphy's BNT toolbox
% function [jtree, cliques, MG, Gtriag] = FindJunctionTree(MG,varargin)
%
% Input:
% MG -> graph
% jt_type = 0 -> block-trees (fastest)
%	    1 -> (default) junction tree using minimum degree elimination order (fast)
%	    2 -> minimum fill elimination order (takes longer)
%
% Output:
% jtree -> tree over cliques
% cliques -> clusters in a junction tree
% MG -> returns the original graph
% Gtriag -> triangulated graph for reference

jt_type = 1;
N = length(MG);
ImpNodes = 1:N;
if ~isempty(varargin)
    jt_type = varargin{1};
    if length(varargin) == 2
        ImpNodes = varargin{2};
    end
    if length(varargin) == 3
        ImpNodes = varargin{2};
        SepSize = varargin{3};
    else
        SepSize = N;
    end
end

MG = MarginalGraph(MG,ImpNodes);
%MG = MG(ImpNodes,ImpNodes);
N = length(MG);

% for constructing block-tree
if (jt_type == 0)
   [jtree, cliques] = construct_block_tree(MG);
   % cliques is a set of NON-OVERLAPPING clusters
   % jtree are edges between the clusters in cliques
   return;
end

MG = setdiag(MG, 0);
Gtriag = MG;
%[S,C] = graphconncomp(G);

stages = {1:N};
ns = ones(size(MG,1),1);
deg = full(sum(MG));

% three cases
% all deg = 0 => empty graph

if (sum(deg) == 0);
    cliques = num2cell(ImpNodes,1);
    jtree = spalloc(N,N,1);
    B = 1;
    return;
end

deg_zero = find(deg == 0);

if isempty(deg_zero)
    if jt_type == 1
        [y,elim_order] = sort(deg,'ascend');
    elseif jt_type == 2
        elim_order = best_first_elim_order(MG, ns, stages);
    end
    [cliques, Gtriag]  = triangulate(MG, elim_order,ImpNodes);
    [jtree, B] = CliquestoJTree(cliques);
else
    % if deg_zero is not empty do the following steps
    cliques = num2cell(deg_zero,1);
    if jt_type == 1
        [y,elim_order] = sort(deg,'ascend');
    elseif jt_type == 2
        elim_order = best_first_elim_order(MG, ns, stages);
    end
    
    [cliques, Gtriag]  = triangulate(MG, elim_order,ImpNodes,cliques);
    
    [jtree, B] = CliquestoJTree(cliques);
end

% find all edges that have a separator of size > SepSize
% and then combine cli        edgesJT = triu(JT.edges); % all edges

edgesJT = triu(jtree); % all edges
indclique = double(edgesJT > SepSize);
[S,C] = graphconncomp(indclique+indclique');

if ~isempty(S)
    for k = 1:S
        cc = find(C == k); % nodes in this component
        new_clique{k} = [];
        for dd = 1:length(cc)
            new_clique{k} = union(new_clique{k},cliques{cc(dd)});
        end
    end
    clear cliques jtree;
    cliques = new_clique;
    [jtree, B] = CliquestoJTree(cliques);
end
end
