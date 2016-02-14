% Conditional Independence Based Cross-Validation

function G = CITCrossValidation(X,G,etaSize)

% Performs conditional independence tests
% to determine if edge (i,j) in G
% Input
% G: graph with edges
% X: nxp observation matrix
%
% Output
% final estimated graph G

G = setdiag(G,0);
p = size(X,2);
if nargin == 2
    eta = max(0,p-2);
else
    eta = etaSize;
end
% if eta is large, only check for all nodes

if eta < 1
    for i = 0:eta
        indedges = find(triu(G,1) == 1); % list of edges to search over
        for e = 1:length(indedges)
            [ii jj] = ind2sub([p p], indedges(e));
            if (i > 0)
                SepSearch = setdiff(1:p,[ii jj]); % list of separators to search over
                if ~isempty(SepSearch)
                    comb = combnk(SepSearch,i); % all possible separators
                else
                    comb = [];
                end
                for j = 1:size(comb,1) % search over all possible separators
                    S = comb(j,:);
                    bool = TestEdgeCV(ii,jj,S,X); % bool = 1 => edge present
                    if (bool == 0)
                        G(ii,jj) = 0; G(jj,ii) = 0;
                        break;
                    end
                end
            else
                bool = TestEdgeCV(ii,jj,[],X);
                if (bool == 0)
                    G(ii,jj) = 0; G(jj,ii) = 0;
                end
            end
        end
    end
else
    indedges = find(triu(G,1) == 1); % list of edges to search over
    for e = 1:length(indedges)
        [ii jj] = ind2sub([p p], indedges(e));
        S = setdiff(1:p,[ii jj]); % list of separators to search over
%         nei = FindNeighborsUndirected(H,i,1);
%         nej = FindNeighborsUndirected(H,j,1);
%         Sij = SeparatorSearchSpace(nei,nej,i,j,H,pcType,k);
        bool = TestEdgeCV(ii,jj,S,X); % bool = 1 => edge present
        if (bool == 0)
            G(ii,jj) = 0; G(jj,ii) = 0;
        end
    end
end

end

function bool = TestEdgeCV(ii,jj,S,X)

K = 5;
n = size(X,1);
numTimes = 100; % number of times 10-fold CV is performed
bool = 0;
if isempty(S)
    for ct = 1:numTimes
        indcv = crossvalind('Kfold',n,K);
        loglike0 = 0;
        loglike1 = 0;
        for k = 1:K
            Xtrain = X(indcv ~= k,[ii jj]);
            Xtest = X(indcv == k,[ii jj]);
            Strain = cov(Xtrain);
            Stest = cov(Xtest);
            loglike0 = loglike0 + -1 * log(prod(diag(Strain))) ...
                - trace(inv([Strain(1,1) 0 ; 0 Strain(2,2)]) * Stest);
            loglike1 = loglike1 + -1 * log(det(Strain)) ...
                - trace(inv(Strain) * Stest);
        end
        if loglike1 > loglike0
            bool = bool + 1;
        else
            break;
        end
    end
else
    for ct = 1:numTimes
        indcv = crossvalind('Kfold',n,K);
        loglike0 = 0;
        loglike1 = 0;
        for k = 1:K
            clear CondTest;
            Xtrain = X(indcv ~= k,[ii jj S]);
            Xtest = X(indcv == k,[ii jj S]);
            Strain = cov(Xtrain);
            Xi = Xtest(:,1);
            Xj = Xtest(:,2);
            Xk = Xtest(:,[3:end]); % n x |S|
            Sijgk = Strain(1:2,1:2) - Strain(1:2,3:end) * inv(Strain(3:end,3:end)) * Strain(1:2,3:end)';
            condMean = Strain(1:2,3:end) * inv(Strain(3:end,3:end)) * Xk'; % 2 x n
            CondTest(:,1) = Xi - condMean(1,:)';
            CondTest(:,2) = Xj - condMean(2,:)';
            loglike1 = loglike1 + -1 * log(det(Sijgk)) - trace(inv(Sijgk) * cov(CondTest));
            
            Sigk = Strain(1,1) - Strain(1,3:end) * inv(Strain(3:end,3:end)) * Strain(1,3:end)';
            condMean = Strain(1,3:end) * inv(Strain(3:end,3:end)) * Xk'; % 1 x n
            CondTest1 = Xi - condMean(1,:)';
            loglike0 = loglike0 + -1 * log(det(Sigk)) - trace(inv(Sigk) * cov(CondTest1));
            
            Sjgk = Strain(2,2) - Strain(2,3:end) * inv(Strain(3:end,3:end)) * Strain(2,3:end)';
            condMean = Strain(2,3:end) * inv(Strain(3:end,3:end)) * Xk'; % 1 x n
            CondTest2 = Xj - condMean(1,:)';
            loglike0 = loglike0 + -1 * log(det(Sjgk)) - trace(inv(Sjgk) * cov(CondTest2));
        end
        if loglike1 > loglike0
            bool = bool + 1;
        else
            break;
        end
    end
end

bool = bool/numTimes > 0.99;

end
