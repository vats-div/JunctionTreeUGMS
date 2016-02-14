function G = UGMS_GLasso(X,H,Hest,options)

% function [G, lambdaChosen] = UGMS_GLasso(X,H,Hest,options)
%
% Only estimate edges over Hest using the graph H that contains a superset
% of the true edges.
%
% options.lambda        -> regularization parameter
% options.EstimateLambda-> if 1, then use model selection else use default
%                           lambda or the lambda given in options.lambda

% options.METHOD        -> Specifies model selection method
%                          BIC(1), STAR(2), ORACLE(3), EBIC(4), EDGES(5)
%                          1 -> BIC
%                          2 -> STAR (proposed by Han Liu et. al.)
%                          3 -> ORACLE, assumes true graph is known
%                          4 -> extended BIC with some specified value of
%                               gamma (see options.gamma)
%                          5 -> Estimate a graph with a certain number of
%                               edges given by options.NumEdges
% options.NumRegul      -> Total number of regularization parameters
% options.SearchLambda  -> Search space for lambda
% options.SearchFurther -> if 1, then search for lambda twice
% options.gamma         -> If using EBIC, then specify gamma
% options.Nodes         -> Nodes over which GLasso should be applied
%                          \bar{R} in the paper
% options.Gtrue         -> true graph to be estimated
%
% Author: Divyanshu Vats, vats.div@gmail.com
%
% Reference
% D. Vats and R. D. Nowak, "A Junction Tree Framework for Undirected
% Graphical Model Selection"
%
% author: Divyanshu Vats


p = size(X,2);
lambda = 0.2;
useMethod = 1;
displayInd = 0;
Hest = setdiag(Hest,0);
numEdgesToEstimate = sum(Hest(:))/2;

% If H is not specified, use complete graph
if isempty(H)
    H = ones(p,p);
end

% If Hest is not specified, estimate all edges  in H
if isempty(Hest);
    Hest = H;
end

H = setdiag(H,0);
Hest = setdiag(Hest,0);

% READ OPTIONS

if isfield(options,'SearchFurther')
    SearchFurther = options.SearchFurther;
else
    SearchFurther = 0;
end

if isfield(options,'displayInd')
    displayInd = options.displayInd;
end

mydisplay(['NumEdges = ' num2str(numEdgesToEstimate)],displayInd);


if isfield(options,'Nodes')
    cc = options.Nodes;
else
    [indi indj] = find(triu(Hest) == 1);
    cc = unique([indi indj]);
end

if isfield(options,'lambda')
    lambda = options.lambda;
else
    mydisplay('No lambda specified, using 0.2 or searching for it',displayInd);
end

% Indicator of whether model selection algorithms need to be used.
if isfield(options,'EstimateLambda')
    EstimateLambda = options.EstimateLambda;
end

% If EstimateLambda = 0, use the given lambda to estimate graph
if EstimateLambda == 0
    G = EstimateGraphgLasso(X,H,Hest,lambda,cc);
    return;
end

% if EstimateLambda = 1, find which model selection algorithm to use
if isfield(options,'METHOD')
    useMethod = options.METHOD;
end

% Number of regularization paramters
if isfield(options,'NumRegul')
    NumRegul = options.NumRegul;
else
    NumRegul = 10;
end

% Search space for regularization parameter
if isfield(options,'SearchLambda')
    SearchLambda = options.SearchLambda;
    NumRegul = length(SearchLambda);
else
    SearchLambda = linspace(1.0,0.01,NumRegul);
end

SearchLambda = sort(SearchLambda,'ascend');

% Read the true graph if given
if isfield(options,'Gtrue')
    Gtrue = options.Gtrue;
else
    Gtrue = zeros(p,p);
    options.Gtrue = Gtrue;
end

% end of reading options

% if there are only a few edges to be estimating, enumerate them all
% if the total number of nodes is small, use hypothesis testing to test for
% edges
if useMethod ~= 3
    if sum(Hest(:))/2 < 7
        useMethod = 6;
    end

    if length(cc) < 8
        useMethod = 7;
    end

end

switch useMethod
    case 1
        mydisplay('Using BIC',displayInd);
        gamma = 0;
        G = FindGraphUsingBIC(X,H,Hest,SearchLambda,cc,SearchFurther,gamma,displayInd);
    case 2
        mydisplay('Using STAR',displayInd);
        if isfield(options,'beta')
            beta = options.beta;
        else
            beta = 0.05;
        end
        G = FindGraphUsingSTAR(X,H,Hest,SearchLambda,cc,beta,SearchFurther,displayInd);
    case 3
        mydisplay('Using Oracle',displayInd);
        G = FindGraphUsingOracle(X,H,Hest,SearchLambda,cc,Gtrue,SearchFurther,displayInd);
    case 4
        mydisplay('Using EBIC',displayInd);
        if isfield(options,'gamma')
            gamma = options.gamma;
        else
            gamma = 0.5;
        end
        G = FindGraphUsingBIC(X,H,Hest,SearchLambda,cc,SearchFurther,gamma,displayInd);
    case 5
        mydisplay('Using Edges',displayInd);
        if isfield(options,'NumEdges')
            NumEdges = options.NumEdges;
        else
            NumEdges = p;
        end
        G = FindGraphUsingEDGES(X,H,Hest,cc,NumEdges);
    case 6
        mydisplay('Using score',displayInd);
        G = FindGraphByEnumerating(X,H,Hest,cc,options.gamma);
    case 7
        mydisplay('Using hypothesis tests',displayInd);
        G = CheckEdgesRobust(X,H,Hest,cc);
    otherwise
end

end


function [G,lambdaChosen] = FindGraphUsingEDGES(X,H,Hest,cc,NumEdges)

% vary lambda until there are NumEdges number of edges in G
CurrentNumEdges = sum(Hest(:))/2;
lambdaChosen = 0;
if CurrentNumEdges <= NumEdges
    G = Hest;
    lambdaChosen = 0;
    return;
end

% Start with small lambda and then increase it until number of edges are
% greater than NumEdges

lambda = 0.1;
G = EstimateGraphgLasso(X,H,Hest,lambda,cc);
CurrentNumEdges = sum(G(:))/2;
while CurrentNumEdges < NumEdges
    lambda = lambda/1.2;
    G = EstimateGraphgLasso(X,H,Hest,lambda,cc);
    CurrentNumEdges = sum(G(:))/2;
    if CurrentNumEdges == NumEdges
        lambdaChosen = 0;
        return;
    end
end

lambdaPrev = lambda;
% start with lambda and increase it until |G| < NumEdges
while CurrentNumEdges > NumEdges
    lambdaPrev = lambda;
    lambda = lambdaPrev * 1.1;
    G = EstimateGraphgLasso(X,H,Hest,lambda,cc);
    CurrentNumEdges = sum(G(:))/2;
    
end

% Now Search between lambdaNext and lambdaPrev
SearchLambda = linspace(lambdaPrev,lambda,20);

for lambda = SearchLambda
    G = EstimateGraphgLasso(X,H,Hest,lambda,cc);
    if sum(G(:))/2 < NumEdges
        lambdaChosen = lambda;
        return;
    end
end

end

% Find graph using stability selection
function [G, lambdaChosen] = FindGraphUsingSTAR(X,H,Hest,SearchLambda,cc,beta,SearchFurther,displayInd)

numRegul = round(length(SearchLambda));
NumEdges = sum(Hest(:))/2;
SearchLambda = sort(SearchLambda,'descend');
n = size(X,1);
p = size(X,2);
indl = 0;
b = floor(10 * sqrt(n));
if b >= n
    b = 0.7 * n;
end
N = 20;
% Gest = zeros(p,p);
% options.EstimateLambda = 0;

DD = zeros(1,length(SearchLambda));
Dbar = DD;

for lambda = SearchLambda
    Gest = zeros(p,p);
    mydisplay(['Checking lambda = ' num2str(lambda)],displayInd);
    for k = 1:N
        tt = randperm(n);
        XX = X(tt(1:b),:);
        Gest = Gest + EstimateGraphgLasso(XX,H,Hest,lambda,cc);
    end
    Gest = Gest/N;
    Q = 2 * Gest .* (1 - Gest);
    indl = indl + 1;
    DD(indl) = sum(sum(triu(Q)))/NumEdges;
    Dbar(indl) = max(DD(1:indl));
end
ss = (Dbar < beta) .* (1:length(Dbar));
ss = setdiff(ss,0);
if isempty(ss)
    indMax = 1;
else
    indMax = max(ss);
end

if SearchFurther == 1
    mydisplay(['Searching Around ' num2str(SearchLambda(indMax))],displayInd);
    % Search around indMax
    setSpace = 0;
    indl = 0;
    if indMax == 1
        SearchLambda = linspace(1.2,SearchLambda(indMax),numRegul);
        setSpace = 1;
    end
    if indMax == length(SearchLambda);
        minss = SearchLambda(end) * 0.25;
        SearchLambda = linspace(SearchLambda(2),minss,numRegul);
        setSpace = 1;
    end
    if setSpace == 0
        SearchLambda = linspace(SearchLambda(indMax-1),SearchLambda(indMax+1),numRegul);
    end
    
    DD = zeros(1,length(SearchLambda));
    Dbar = DD;
    
    for lambda = SearchLambda
        Gest = zeros(p,p);
        mydisplay(['Checking lambda = ' num2str(lambda)],displayInd);
        for k = 1:N
            tt = randperm(n);
            XX = X(tt(1:b),:);
            Gest = Gest + EstimateGraphgLasso(XX,H,Hest,lambda,cc);
        end
        Gest = Gest/N;
        Q = 2 * Gest .* (1 - Gest);
        indl = indl + 1;
        DD(indl) = sum(sum(triu(Q)))/NumEdges;
        Dbar(indl) = max(DD(1:indl));
    end
    ss = (Dbar < beta) .* (1:length(Dbar));
    ss = setdiff(ss,0);
    if isempty(ss)
        indMax = 1;
    else
        indMax = max(ss);
    end
end

mydisplay(['Choosing lambda = ' num2str(SearchLambda(indMax))],displayInd);
G = EstimateGraphgLasso(X,H,Hest,SearchLambda(indMax),cc);
lambdaChosen = SearchLambda(indMax);

end

function [G, lambdaChosen] = FindGraphUsingBIC(X,H,Hest,SearchLambda,cc,SearchFurther,gamma,displayInd)

indl = 0;
numRegul = floor(length(SearchLambda));

for lambda = SearchLambda
    indl = indl + 1;
    [GraphEst{indl}, ScoreBIC(indl)] = EstimateGraphgLassoBIC(X,H,Hest,lambda,cc,gamma); %#ok<*AGROW>
end

% find the right lambda
indMin = ceil(argmin(ScoreBIC));

% Further search around SearchLambda(indMin)

if SearchFurther == 1
    mydisplay(['Searching Around ' num2str(SearchLambda(indMin))],displayInd);
    clear ScoreBIC GraphEst;
    % Search around indMin
    setSpace = 0;
    if indMin == 1
        minss = SearchLambda(1) * 0.25;
        SearchLambda = linspace(minss,SearchLambda(2),numRegul);
        setSpace = 1;
    end
    if indMin == length(SearchLambda);
        SearchLambda = linspace(SearchLambda(indMin),1.5,numRegul);
        setSpace = 1;
    end
    if setSpace == 0
        SearchLambda = linspace(SearchLambda(indMin-1),SearchLambda(indMin+1),numRegul);
    end
    
    indl = 0;
    for lambda = SearchLambda
        indl = indl + 1;
        [GraphEst{indl}, ScoreBIC(indl)] = EstimateGraphgLassoBIC(X,H,Hest,lambda,cc,gamma);
    end
end

indMin = ceil(argmin(ScoreBIC));
mydisplay(['Choosing lambda = ' num2str(SearchLambda(indMin))],displayInd);
G = GraphEst{indMin};
lambdaChosen = SearchLambda(indMin);
clear GraphEst;

end

function [G, lambdaChosen] = FindGraphUsingOracle(X,H,Hest,SearchLambda,cc,Gtrue,SearchFurther,displayInd)

indl = 0;
numRegul = floor(length(SearchLambda)*4);

for lambda = SearchLambda
    mydisplay(['Checking lambda = ' num2str(lambda)],displayInd);
    indl = indl + 1;
    GraphEst{indl} = EstimateGraphgLasso(X,H,Hest,lambda,cc);
    GTemp = GraphEst{indl};
%    HTemp = H(cc,cc);
%    HTemp = HTemp - Hest(cc,cc) + GTemp(cc,cc);
    ScoreOracle(indl) = CompareGraphs(GTemp(cc,cc),Gtrue(cc,cc).*Hest(cc,cc));    
end

indMin = argmin(ScoreOracle);

if SearchFurther == 1
    clear ScoreOracle GraphEst;
    
    % Search around indMin
    setSpace = 0;
    if indMin == 1
        minss = SearchLambda(1) * 0.25;
        SearchLambda = linspace(minss,SearchLambda(2),numRegul);
        setSpace = 1;
    end
    if indMin == length(SearchLambda)
        SearchLambda = linspace(SearchLambda(indMin),1.5,numRegul);
        setSpace = 1;
    end
    if setSpace == 0
        SearchLambda = linspace(SearchLambda(indMin-1),SearchLambda(indMin+1),numRegul);
    end

    indl = 0;
    for lambda = SearchLambda
        indl = indl + 1;
        GraphEst{indl} = EstimateGraphgLasso(X,H,Hest,lambda,cc);
        GTemp = GraphEst{indl};
        ScoreOracle(indl) = CompareGraphs(GTemp(cc,cc),Gtrue(cc,cc).*Hest(cc,cc));    
        
    end

end

indMin = argmin(ScoreOracle);
G = GraphEst{indMin};
lambdaChosen = SearchLambda(indMin);
end

% Applies graphical Lasso for a particular lambda
function G = EstimateGraphgLasso(X,H,Hest,lambda,cc)

% only returns edges in Hest that are 1

p = size(H,1);
Hcc = MarginalGraph(H,cc);
Domain = ones(p,p);
Domain(Hest == 1) = 0;
Domain = Domain(cc,cc);
Domain(Hcc == 0) = 2;
G = spalloc(p,p,p);

% this will estimate all edges over cc that are in Hest
invCov = covsel(X(:,cc),lambda,1,1,Domain);

G(cc,cc) = abs(invCov) > 1e-4;
G(Hest ~= 1) = 0;
G = sparse(G);

end

function [G, ScoreBIC] = EstimateGraphgLassoBIC(X,H,Hest,lambda,cc,gamma)

p = size(H,1);
Hcc = MarginalGraph(H,cc);
Domain = ones(p,p);
Domain(Hest == 1) = 0;
Domain = Domain(cc,cc);
Domain(Hcc == 0) = 2;
G = spalloc(p,p,p);

% this will estimate all edges over cc that are in Hest
invCov = covsel(X(:,cc),lambda,1,1,Domain);

G(cc,cc) = abs(invCov) > 1e-4;
G(Hest ~= 1) = 0;
G = sparse(G);
HTemp = Hcc - Hest(cc,cc) + G(cc,cc);
ScoreBIC = BICScore(X(:,cc),HTemp,Hest(cc,cc),gamma);

end