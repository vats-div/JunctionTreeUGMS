function [coef,pval] = PartialCorrCoef(CovMat,i,j,S,n)

% compute partial correlation of i and j given S


ssij = CovMat(i,j) - CovMat(i,S) * inv(CovMat(S,S)) * CovMat(S,j);
ssjj = CovMat(j,j) - CovMat(j,S) * inv(CovMat(S,S)) * CovMat(S,j);
ssii = CovMat(i,i) - CovMat(i,S) * inv(CovMat(S,S)) * CovMat(S,i);

coef = ssij/sqrt(ssii*ssjj);

if nargout > 1
    df = max(n - length(S) - 2,0); % this is a matrix for 'pairwise'
    t = sign(coef) .* Inf;
    k = (abs(coef) < 1);
    t(k) = coef(k) ./ sqrt(1-coef(k).^2);
    t = sqrt(df).*t;
    pval = 2*tcdf(-abs(t),df);
end

end