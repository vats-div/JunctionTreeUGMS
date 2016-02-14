function Score = CVScore(X,G,H,inddc)

Score = 0;
G = setdiag(G,1);

% to K-fold cross-validation with
K = 10;
n = size(X,1);

H = triu(H,1);


for tt = 1:length(inddc)
    indd = inddc{tt};
    for k = 1:K
        indTrain = (indd ~= k);
        indTest = (indd == k);
        XTrain = X(indTrain,:);
        XTest = X(indTest,:);
        if sum(sum(G)) == sum(sum(ones(size(X,2))))
            InvSigma = inv(cov(XTrain));
        else
            InvSigma = InverseCovariance(XTrain,G);
        end
        
        % for each edge estimating, compute $
        
        
        Score = Score + log(det(InvSigma)) - ...
            trace(cov(XTest) * InvSigma);
    end
end



end