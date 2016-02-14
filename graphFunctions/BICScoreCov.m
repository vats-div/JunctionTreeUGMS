function BIC = BICScore(X,H,Hest,InvSigma,gammaGiven)

% H : Estimate of the graph
% Hest : Edges in H that were estimated
% Thus, if Hest(i,j) = 1, then the edge in H(i,j) is estimated using X

% NumParam = sum(H(:).*Hest(:))/2; 
% BIC score = -2 * loglikelihood + |E| log(n) + 4 |E| gamma log(p)

if nargin == 4
	gammaGiven = 0.5;
end
H = setdiag(H,0);
H = H + eye(size(H));
[indi indj] = find(triu(Hest) == 1);


k = sum(H(:).*Hest(:))/2; % total number of edges estimates from Hest
n = size(X,1);
p = size(X,2);

Nodes = (unique([indi indj]));

gamma = gammaGiven;


% Find the inverse covariance matrix over all cc
S = cov(X);

% Likelihood
logL = 0.5 * n * (log(det(InvSigma)) - trace(S * InvSigma));
BIC = -2*logL + k*log(n) + 4*k*gamma*log(length(Nodes));

% NodesC = setdiff(1:p,Nodes);
% InvSigmaNodes = InvSigma(Nodes,Nodes) - InvSigma(Nodes,NodesC) * inv(InvSigma(NodesC,NodesC)) * InvSigma(NodesC,Nodes);
% SNodes = S(Nodes,Nodes) - S(Nodes,NodesC) * inv(S(NodesC,NodesC)) * S(NodesC,Nodes);
% logL = 0.5 * n * (log(det(InvSigmaNodes)) - trace(SNodes * InvSigmaNodes));
% BIC = -2*logL + k*log(n) + 4*k*gamma*log(length(Nodes));


end
