function [ Tree,Cost ] =  UndirectedMaximumSpanningTree (CostMatrix,varargin)
% The function takes CostMatrix as input and returns the maximum spanning tree T
% Uses Kruskal's Algorithm
% Extract the edge weights from the cost matrix
% Sort the edges in a non decreasing order of weights 
% This algorithm is revised by Lowell Guangdi at 2009/06/11.
% revised by Divyanshu Vats on Dec. 8, 2011

n = size (CostMatrix,1); %Number of vertices

if (isempty(varargin))
    s = n - 1;
else
    s = min(varargin{1},n-1);
    if (s == 0)
        Tree = eye(n,n); Cost = 0;
        return;
    end
end


EdgeWeights = 0;         %Edges and corresponding weights
EdgeWeightsCounter = 0;

ind_edges = find(triu(CostMatrix,1) > 0);
[ii, jj] = ind2sub(size(CostMatrix),ind_edges);
EdgeWeights = zeros(length(ind_edges),3);
for d = 1:length(ind_edges)
    i = ii(d); j = jj(d);
    EdgeWeights(d,1) = CostMatrix(i,j);
    EdgeWeights(d,2) = i;
    EdgeWeights(d,3) = j; 
end

if (EdgeWeights == 0)
    Tree = eye(n,n); Cost = 0;
    return;
end

SortedEdgeWeights = 0;
SortedEdgeWeights = sortrows(EdgeWeights);
% First column of SortedEdgeWeights are the weights
% Second and third column are the vertices that the edges connect
m = size(SortedEdgeWeights,1); % number of edges 

% We use the Disjoint sets data structures to detect cycle while adding new
% edges. Union by Rank with path compression is implemented here.

% Assign parent pointers to each vertex. Initially each vertex points to 
% itself. Now we have a conceptual forest of n trees representing n disjoint 
% sets 
global ParentPointer ;
ParentPointer = 0;
ParentPointer(1:n) = 1:n;

% Assign a rank to each vertex (root of each tree). Initially all vertices 
% have the rank zero.
TreeRank = 0;
TreeRank(1:n) = 0;

% Visit each edge in the sorted edges array
% If the two end vertices of the edge are in different sets (no cycle), add
% the edge to the set of edges in maximum spanning tree
MSTreeEdges = 0;
MSTreeEdgesCounter = 0; i = m;
while ((MSTreeEdgesCounter < s) && (i>=1))
%Find the roots of the trees that the selected edge's two vertices
%belong to. Also perform path compression.
    root1=0; root2=0; temproot=0;
    temproot = SortedEdgeWeights(i,2);
    root1 = FIND_PathCompression(temproot);
  
    temproot = SortedEdgeWeights(i,3);
    root2 = FIND_PathCompression(temproot);
    
    if (root1 ~= root2)
        MSTreeEdgesCounter = MSTreeEdgesCounter + 1;
        MSTreeEdges(MSTreeEdgesCounter,1:3) = SortedEdgeWeights(i,:);
        if (TreeRank(root1)>TreeRank(root2))
            ParentPointer(root2)=root1;
        else
            if (TreeRank(root1)==TreeRank(root2))
               TreeRank(root2)=TreeRank(root2) + 1;
            end
            ParentPointer(root1)=root2;
        end
    end
    i = i - 1;
end
Num_Edges = MSTreeEdgesCounter;
MSTreeEdgesCounter = 0;
Tree = 0;
Tree(1:n,1:n)=0;
while (MSTreeEdgesCounter < Num_Edges)
    MSTreeEdgesCounter = MSTreeEdgesCounter + 1;
    Tree(MSTreeEdges(MSTreeEdgesCounter,2),MSTreeEdges(MSTreeEdgesCounter,3))=1;
    Tree(MSTreeEdges(MSTreeEdgesCounter,3),MSTreeEdges(MSTreeEdgesCounter,2))=1;
end
%T

Cost = 0;


end

function [parent] = FIND_PathCompression(temproot)

global ParentPointer;
ParentPointer(temproot);
if (ParentPointer(temproot)~=temproot)
    ParentPointer(temproot) = FIND_PathCompression(ParentPointer(temproot));
end
parent = ParentPointer(temproot);
end