function G = HubGraph(p)

% choose components of size between
% 2 and sqrt(p)
psz = ceil(sqrt(p));
compsize = 1 + ceil( psz * rand(p,1));
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

rho = 1/(compsize(1)+1);
G = zeros(compsize(1),compsize(1));
G(1,:) = rho; G(:,1) = rho;

for k = 2:length(compsize)
    rho = 1/(compsize(k)+1);
    Gb = zeros(compsize(k),compsize(k));
    Gb(1,:) = rho; Gb(:,1) = rho;
    G = blkdiag(G,Gb);
end

G = setdiag(G,1);

end