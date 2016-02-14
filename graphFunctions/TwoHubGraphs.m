function [invSigma, Gtrue] = TwoHubGraphs(p,d1,d2,perc)

% smaller part has hubs of size d1
% larger part has hubs of size d2

p1 = perc * p; p2 = p - p1;
Gtrue1 = myHubGraph(p1,d1);
Gtrue2 = myHubGraph(p2,d2);

invSigma = zeros(p,p);
invSigma(1:p1,1:p1) = Gtrue1;
invSigma(p1+1:end,p1+1:end) = Gtrue2;
Gtrue = double(abs(invSigma) > 0);
Gtrue = setdiag(Gtrue,0);

end

function G = myHubGraph(p,d)

% choose components of size d

compsize = d * ones(p,1);
ss = (cumsum(compsize) <= p) .* (1:p)';
indmax = max(ss);

if sum(compsize(1:indmax)) < p
    compsize(indmax+2:end) = [];
    compsize(indmax+1) = p - sum(compsize(1:indmax));
end

if sum(compsize(1:indmax)) > p
    compsize(indmax:end) = [];
    compsize(indmax) = p - sum(compsize(1:indmax-1));
end

if sum(compsize(1:indmax)) == p
    compsize(indmax+1:end) = [];
end

G = zeros(compsize(1),compsize(1));
rho = 1/(d-1) - 0.0000001;

G(1,:) = rho; G(:,1) = rho;

for k = 2:length(compsize)
    Gb = zeros(compsize(k),compsize(k));
    Gb(1,:) = rho; Gb(:,1) = rho;
    G = blkdiag(G,Gb);
end

G = setdiag(G,1);

end