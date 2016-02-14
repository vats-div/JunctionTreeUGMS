function G = UGMS_PC(X,H,Hest,options) %#ok<*INUSD,*INUSL>

% function [G SearchComp] = UGMS-PC(Delta,X,H,Hest,options)
%
% Implements the PC-Algorithm for Undirected Graphical Model Selection when
% the graphical model is Gaussian
% Can be extended to discrete by changing the conditional independence test
% H: graph that contains the true edges (or most of the true edges)
% Hest: Only estimate edges in Hest
%
% options.Delta     -> Separator Search Space (should be a vector,like 0:2)
% options.xi        -> threshold for conditional independence
% options.Nodes     -> Nodes in H over which edges are estimated
% options.pcType    -> 1 to use marginal graphs and 2 to use neighbors
% options.EstimateXi-> 1 Use model selection to find xi
% options.METHOD    -> BIC(1) or STAR(2) or ORACLE(3) or EBIC(4), default is BIC
% options.NumRegul  -> Number of Parameters to use
% options.SearchXi  -> Search Space

% Initialization
p = size(X,2); % number of nodes;
n = size(X,1);
xi = 0.2;
Delta = 0:1;
pcType = 1;
EstimateXi = 0;
useMethod = 1;
xiChosen = 0;
displayInd = 0;

numEdgesToEstimate = sum(Hest(:))/2;

% Preprocessing
if isempty(H)
    H = ones(p,p);
end

if isempty(Hest)
    Hest = H;
end

% Make diagonal zero
Hest = setdiag(Hest,0);
H = setdiag(H,0);

% Get options
if isfield(options,'Delta')
        Delta = options.Delta;
end

if isfield(options,'displayInd')
    displayInd = options.displayInd;
end

mydisplay(['NumEdges = ' num2str(numEdgesToEstimate)],displayInd);


if isfield(options,'SearchFurther')
        SearchFurther = options.SearchFurther;
else
        SearchFurther = 1;
end

if isfield(options,'xi')
    xi = options.xi;
end

if isfield(options,'Nodes')
    Nodes = options.Nodes; %#ok<*NASGU>
else
    Nodes = 1:p;
    options.Nodes = Nodes;
end

if isfield(options,'pcType')
    pcType = options.pcType;
end

if isfield(options,'EstimateXi')
    EstimateXi = options.EstimateXi;
end
% end of options


% Initialization
G = Hest;

% if no edges in G, do nothing
if all(G(:))
    return;
end

% Do not estimate Xi, use given value of xi
if EstimateXi == 0
    G = MainUGMSPCMain(Delta,X,H,G,xi,pcType,options);
    G = sparse(G);
    return;
end

% If EstimateXi == 1, then search for xi using model selection

if isfield(options,'METHOD')
    useMethod = options.METHOD;
end

if isfield(options,'gamma')
    gamma = options.gamma;
else
    gamma = 0.5;
end

% Number of Parameters to choose
if isfield(options,'NumRegul')
    NumRegul = options.NumRegul;
else
    NumRegul = 10;
end

if isfield(options,'SearchXi')
    SearchXi = options.SearchXi;
    NumRegul = length(SearchXi);
else
    SearchXi = linspace(0.5,0.1,NumRegul);
end

indxi = 0;

%display('Starting Search of a threshold...');

if isfield(options,'Gtrue')
    Gtrue = options.Gtrue;
else
    Gtrue = zeros(p,p);
    options.Gtrue = Gtrue;
    
end

SearchXi = sort(SearchXi,'ascend');

% if there are only a few edges to be estimating, enumerate them all
% if the total number of nodes is small, use hypothesis testing to test for
% edges
%if useMethod ~= 3
    if sum(Hest(:))/2 < 7
        useMethod = 6;
    end

    if length(Nodes) < 8
        useMethod = 7;
    end

%end

switch useMethod
    case 1
        mydisplay('Using BIC',displayInd);
        gamma = 0;
        G = FindGraphUsingBIC(Delta,X,H,Hest,SearchXi,pcType,Nodes,SearchFurther,gamma);
    case 2
        mydisplay('Using StaR',displayInd);
        G = FindGraphUsingStaR(Delta,X,H,Hest,options,SearchXi);
    case 3
        mydisplay('Using Oracle',displayInd)
        G = FindGraphUsingOracle(Delta,X,H,Hest,SearchXi,pcType,Nodes,Gtrue,SearchFurther);
    case 4
        mydisplay('Using EBIC',displayInd);
        if isfield(options,'gamma')
            gamma = options.gamma;
        else
            gamma = 0.5;
        end
        G = FindGraphUsingBIC(Delta,X,H,Hest,SearchXi,pcType,Nodes,SearchFurther,gamma);
    case 5
        mydisplay('Match NumEdges',displayInd);
        if isfield(options,'NumEdges')
            NumEdges = options.NumEdges;
        else
            NumEdges = p;
        end
        G = FindGraphUsingEDGES(X,H,Hest,NumEdges,pcType,options);
        
    case 6
        mydisplay('Using score',displayInd);
        xiChosen = NaN;
        G = FindGraphByEnumerating(X,H,Hest,Nodes,options.gamma);
    case 7
        mydisplay('HT',displayInd);
        xiChosen = NaN;
        G = CheckEdgesRobust(X,H,Hest,Nodes);
    otherwise
        G = FindGraphUsingBIC(Delta,X,H,Hest,options,SearchFurther);
end

end

function G = FindGraphUsingEDGES(X,H,Hest,NumEdges,pcType,options)

Delta = options.Delta;

CurrentNumEdges = sum(Hest(:))/2;
if CurrentNumEdges <= NumEdges
    G = Hest;
    return;
end

% Start with small lambda and then increase it until number of edges are
% greater than NumEdges

lambda = 0.2;
G = MainUGMSPC(Delta,X,H,Hest,lambda,pcType,options);
CurrentNumEdges = sum(G(:))/2;
while CurrentNumEdges < NumEdges
    lambda = lambda/1.5;
    G = MainUGMSPC(Delta,X,H,Hest,lambda,pcType,options);
    CurrentNumEdges = sum(G(:))/2;
    if CurrentNumEdges == NumEdges
        return;
    end
end

lambdaPrev = lambda;
% start with lambda and increase it until |G| < NumEdges
while CurrentNumEdges > NumEdges
    lambdaPrev = lambda;
    lambda = lambdaPrev * 1.5;
    G = MainUGMSPC(Delta,X,H,Hest,lambda,pcType,options);
    CurrentNumEdges = sum(G(:))/2;
end

% Now Search between lambdaNext and lambdaPrev
SearchLambda = linspace(lambdaPrev,lambda,10);

for lambda = SearchLambda
    G = MainUGMSPC(Delta,X,H,Hest,lambda,pcType,options);
    if sum(G(:))/2 < NumEdges
        lambdaChosen = lambda;
        return;
    end
end

end


function G = FindGraphUsingBIC(Delta,X,H,Hest,SearchXi,pcType,Nodes,SearchFurther,gamma)

indxi = 0;
numRegul = round(length(SearchXi));
options.Nodes = Nodes;
H = double(H > 0);

% Use BIC
for xi = SearchXi
    % First check for independence
    %display(['Checking xi = ' num2str(xi)]);
    indxi = indxi + 1;
    GraphEst{indxi} = MainUGMSPC(Delta,X,H,Hest,xi,pcType,options); % estimates edges over Hest
    GTemp = GraphEst{indxi};
    HTemp = MarginalGraph(H,Nodes);
    HTemp = HTemp - Hest(Nodes,Nodes) + GTemp(Nodes,Nodes);
    ScoreBIC(indxi) = BICScore(X(:,Nodes),HTemp,Hest(Nodes,Nodes),gamma);
    %if sum(sum(GraphEst{indxi})) == 0
    %	break;
    %    end
end

% choose graph with minimum score

indMin = argmin(ScoreBIC);
if (SearchFurther == 1 && sum(sum(GraphEst{indMin})) > 0)
    clear ScoreBIC GraphEst;
    %display(['Searching Around ' num2str(SearchXi(indMin))]);
    setSpace = 0;
    if indMin == 1
        minsearch = SearchXi(1) * 0.5;
        SearchXi = linspace(minsearch,SearchXi(2),numRegul);
        setSpace = 1;
    end
    if indMin == length(SearchXi)
        SearchXi = linspace(SearchXi(indMin),1.0,numRegul);
        setSpace = 1;
    end
    if setSpace == 0
        SearchXi = linspace(SearchXi(indMin-1),SearchXi(indMin+1),numRegul);
    end

    indxi = 0;
    for xi = SearchXi
        % First check for independence
        %display(['Checking xi = ' num2str(xi)]);
        indxi = indxi + 1;
        GraphEst{indxi} = MainUGMSPC(Delta,X,H,Hest,xi,pcType,options); % will estimate edges over Hest
    end

    for k = 1:length(SearchXi)
        GTemp = GraphEst{k};
        %HTemp = H(Nodes,Nodes);
        HTemp = MarginalGraph(H,Nodes);
        HTemp = HTemp - Hest(Nodes,Nodes) + GTemp(Nodes,Nodes);
        ScoreBIC(k) = BICScore(X(:,Nodes),HTemp,Hest(Nodes,Nodes));
    end
    indMin = argmin(ScoreBIC);
end

%display(['Choosing xi = ' num2str(SearchXi(indMin))]);
G = GraphEst{indMin};

end

function G = FindGraphUsingOracle(Delta,X,H,Hest,SearchXi,pcType,Nodes,Gtrue,SearchFurther)

%display('here');
indxi = 0;
p = size(X,2);
numRegul = round(length(SearchXi));

options.Nodes = Nodes;

% Use BIC
for xi = SearchXi
    % First check for independence
    %display(['Checking xi = ' num2str(xi)]);
    indxi = indxi + 1;
    GraphEst{indxi} = MainUGMSPC(Delta,X,H,Hest,xi,pcType,options); % will estimate edges over Hest
end

H = double(H > 0);

for k = 1:length(SearchXi)
    GTemp = GraphEst{k};
    HTemp = H(Nodes,Nodes);
    HTemp = HTemp - Hest(Nodes,Nodes) + GTemp(Nodes,Nodes);
    ScoreOracle(k) = CompareGraphs(GTemp(Nodes,Nodes),Gtrue(Nodes,Nodes).*Hest(Nodes,Nodes));
end

% choose graph with minimum score
indMin = argmin(ScoreOracle);
if SearchFurther == 1
    clear ScoreOracle GraphEst;
    %display(['Searching Around ' num2str(SearchXi(indMin))]);
    setSpace = 0;
    if indMin == 1
        minsearch = SearchXi(1) * 0.5;
        SearchXi = linspace(minsearch,SearchXi(2),numRegul);
        setSpace = 1;
    end
    
    if indMin == length(SearchXi)
        SearchXi = linspace(SearchXi(indMin),1.0,numRegul);
        setSpace = 1;
    end
    if setSpace == 0
        SearchXi = linspace(SearchXi(indMin-1),SearchXi(indMin+1),numRegul);
    end
    indxi = 0;
    for xi = SearchXi
        % First check for independence
        %display(['Checking xi = ' num2str(xi)]);
        indxi = indxi + 1;
        GraphEst{indxi} = MainUGMSPC(Delta,X,H,Hest,xi,pcType,options); % will estimate edges over Hest
    end
    
    for k = 1:length(SearchXi)
        GTemp = GraphEst{k};
        HTemp = H(Nodes,Nodes);
        HTemp = HTemp - Hest(Nodes,Nodes) + GTemp(Nodes,Nodes);
        ScoreOracle(k) = CompareGraphs(GTemp(Nodes,Nodes),Gtrue(Nodes,Nodes).*Hest(Nodes,Nodes));
    end
    indMin = argmin(ScoreOracle);
end
%display(['Choosing xi = ' num2str(SearchXi(indMin))]);
G = GraphEst{indMin};
end

function G = FindGraphUsingStaR(Delta,X,H,Hest,options,SearchXi)

n = size(X,1);
p = size(X,2);

b = floor(10 * sqrt(n)); % number of samples

if b >= n
    b = floor(0.8 * n);
end

N = 50; % number of times subsample

maxEta = max(Delta);
minEta = min(Delta);

indxi = 0;
if isfield(options,'beta')
    beta = options.beta;
else
    beta = 0.05;
end

N = 30;
Gest = zeros(p,p);
options.EstimateXi = 0;


% Compute all subsamples
for k = 1:N
    tt = randperm(n);
    XX{k} = X(tt(1:b),:);
end

options.EstimateXi = 0;

SearchXi = sort(SearchXi,'descend');

%for eta = Delta
%    options.Delta = 0:eta;
    indxi = 0;
    clear Dbar;
    clear D;
    NumEdges = sum(Hest(:))/2;
    for xi = SearchXi
        options.EstimateXi = 0;
        Gest = zeros(p,p);
        %display(['Checking xi = ' num2str(xi)]);
        options.xi = xi;
        for k = 1:N
            Gest = Gest + UGMS_PC(XX{k},H,Hest,options);
        end
        Gest = Gest/N;
        Q = 2 * Gest .* (1 - Gest);
        indxi = indxi + 1;
        D(indxi) = sum(sum(triu(Q)))/NumEdges;
        Dbar(indxi) = max(D(1:indxi))
        if (Dbar(end) >= beta)
           break;
        end
    end
    
    ss = (Dbar < beta) .* (1:length(Dbar));
    ss = setdiff(ss,0);
    
    if isempty(ss)
        indMax = length(SearchXi);
    else
        indMax = max(ss);
    end
    
    % Search around indMax
    display(['Choosing xi = ' num2str(SearchXi(indMax))]);
    options.xi = SearchXi(indMax);
    G = UGMS_PC(X,H,Hest,options);
    G = double(G);
    G = (abs(G) > 0);
    
    % Update H and Hest
    H = H - Hest + G; % current estimate of H
    Hest = G; % current estimate of edges to estimate in H
%end

end

function G = MainUGMSPC(Delta,X,H,G,xi,pcType,options)

G = MainUGMSPCMain(Delta,X,H,G,xi,pcType,options);

end

function G = MainUGMSPCMain(Delta,X,H,G,xi,pcType,options)

%Nodes = options.Nodes;
p = size(X,2);

if min(Delta) == 0
    cc = mycorr(X);
    bool = (abs(cc) < xi);
    G(bool == 1) = 0;
    H(bool == 1) = 0;
    Delta = setdiff(Delta,0);
end

for k = Delta
    % find all edges that need to estimated
    ind_edges = find(triu(G) == 1);
    %ind_edges = ind_edges(randperm(length(ind_edges)));
    num_edges = length(ind_edges);
    for e = 1:num_edges
        [i j] = ind2sub([p p], ind_edges(e)); % the edge
        nei = FindNeighborsUndirected(H,i,1);
        nej = FindNeighborsUndirected(H,j,1);
        Sij = SeparatorSearchSpace(nei,nej,i,j,H,pcType,k);
        if ~isempty(Sij)
            comb = combnk(Sij,k);
            for ct = 1:size(comb,1)
                S = comb(ct,:);
                % conditional independence
                % Can replace this by PartialCorrCoef(cov(X),i,j,S,n) too
                cc = partialcorr(X(:,i),X(:,j),X(:,S));
                bool = (abs(cc) < xi); %bool = 1 => X_i indep X_j | X_S
                if (bool == 1)
                    G(i,j) = 0; G(j,i) = 0; H(i,j) = 0; H(j,i) = 0;
                    break;
                end
            end
        else
            H(i,j) = 0; H(j,i) = 0;
        end
    end
end

G = setdiag(G,0);

end