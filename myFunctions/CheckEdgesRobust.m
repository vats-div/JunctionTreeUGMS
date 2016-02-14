% Conditional Independence Based Cross-Validation

function Gcc = CheckEdgesRobust(X,G,Hest,cc)

% Test edges in G(cc,cc) independently

p = size(X,2);
Gcc = spalloc(p,p,p);
Gcc(cc,cc) = Hest(cc,cc);
indedges = find(triu(Gcc,1) == 1);

for e = 1:length(indedges)
    [ii jj] = ind2sub([p p], indedges(e));
    S1 = setdiff(cc,[ii jj]);
    bool = TestEdgeCV(ii,jj,S1,X); % bool = 1 => edge present
    if (bool == 0)
        Gcc(ii,jj) = 0; Gcc(jj,ii) = 0;
    end
    if (bool ~= 0)
        S2 = setdiff(FindNeighborsUndirected(G,ii,1),jj);
        bool = TestEdgeCV(ii,jj,S2,X); % bool = 1 => edge present
        if (bool == 0)
            Gcc(ii,jj) = 0; Gcc(jj,ii) = 0;
        end
    end
    if (bool ~= 0)
        S2 = setdiff(FindNeighborsUndirected(G,jj,1),ii);
        bool = TestEdgeCV(ii,jj,S2,X); % bool = 1 => edge present
        if (bool == 0)
            Gcc(ii,jj) = 0; Gcc(jj,ii) = 0;
        end
    end
end


end

function bool = TestEdgeCV(ii,jj,S,X)

K = 5;
n = size(X,1);
numTimes = 50; % number of times 10-fold CV is performed
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
