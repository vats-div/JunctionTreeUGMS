function invSigma = InverseCovariance(X,G)

% function invSigma = InverseCovariance(X,G)
% Finds the inverse of the covariance on a graph G given data X

p = size(X,2);
n = size(X,1);
G = sparse(G);
[S,C] = graphconncomp(G);

if (S == 1)
    GG = triu(G,1) + tril(ones(p,p));
    invSigma = InverseCovarianceAllNodes(X,GG);
else
    invSigma = spalloc(p,p,p);
    for k = 1:S
        cc = find(C == k);
        pc = length(cc);
        Gc = G(cc,cc);
        GG = triu(Gc,1) + tril(ones(pc,pc));
        if (~all(GG(:)))
            invSc = InverseCovarianceAllNodes(X(:,cc),GG);
        else
            invSc = inv(X(:,cc)' * X(:,cc)/n);
        end
        invSigma(cc,cc) = invSc;
    end
end


end

function invSigma = InverseCovarianceAllNodes(X,G)

n = size(X,1);
OPTIONS.sigma      = 1;
OPTIONS.scale_data = 0;
OPTIONS.plotyes = 0;
OPTIONS.tol     = 1e-6;
OPTIONS.smoothing = 1;

S = cov(X);
[Iomega,Jomega] = find(G == 0);
[blk,At,C,b,mu] = genSDPdata(Iomega,Jomega,S,0,0);
[obj,SigmaEst,y,Z,info,runhist] = logdetPPA(blk,At,C,b,mu,OPTIONS);
invSigma = SigmaEst{1};

end
