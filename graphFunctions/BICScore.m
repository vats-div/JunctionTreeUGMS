function BIC = BICScore(X,H,Hest,gammaGiven,InvSigma)

% H : Estimate of the graph
% Hest : Edges in H that were estimated
% Thus, if Hest(i,j) = 1, then the edge in H(i,j) is estimated using X

% NumParam = sum(H(:).*Hest(:))/2;
% BIC score = -2 * loglikelihood + |E| log(n) + 4 |E| gamma log(p)

if nargin == 3
    gammaGiven = 0.5;
end

H = setdiag(H,0);
Hest = setdiag(Hest,0);

[indi indj] = find(triu(Hest) == 1);
NumEdgesToEstimate = sum(Hest(:)) + length(unique([indi indj]));
k = sum(H(:).*Hest(:))/2; % total number of edges estimated from Hest
n = size(X,1);

gamma = gammaGiven;

S = cov(X);
H = H + eye(size(H));
if nargin < 5
    % Find the inverse covariance matrix over all cc
    if sum(sum(H)) == sum(sum(ones(size(X,2))))
        InvSigma = inv(S);
    else
        InvSigma = InverseCovariance(X,H);
    end
end

% Likelihood
logL = 0.5 * n * (log(det(InvSigma)) - trace(S * InvSigma));
if isinf(logL)
    logL = -Inf;
end
%BIC = -2*logL + k*log(n) + 2*k*gamma*log(NumEdgesToEstimate);
BIC = -2*logL + k*log(n) + 4*k*gamma*log(size(Hest,1));

% NodesC = setdiff(1:p,Nodes);
% InvSigmaNodes = InvSigma(Nodes,Nodes) - InvSigma(Nodes,NodesC) * inv(InvSigma(NodesC,NodesC)) * InvSigma(NodesC,Nodes);
% SNodes = S(Nodes,Nodes) - S(Nodes,NodesC) * inv(S(NodesC,NodesC)) * S(NodesC,Nodes);
% logL = 0.5 * n * (log(det(InvSigmaNodes)) - trace(SNodes * InvSigmaNodes));
% BIC = -2*logL + k*log(n) + 4*k*gamma*log(length(Nodes));


end
