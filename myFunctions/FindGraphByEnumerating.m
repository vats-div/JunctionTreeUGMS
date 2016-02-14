function G = FindGraphByEnumerating(X,H,Hest,cc,gamma)

Hm = MarginalGraph(H,cc);
Hest = triu(Hest,1);
HestSmall = Hest(cc,cc);
EdgesToEstimate = find(HestSmall == 1);
HestSmall = HestSmall + HestSmall';
NumEst = length(EdgesToEstimate);

for k = 0:2^NumEst-1
    patt = BinaryVector(k,NumEst);
    HTemp = triu(Hm,1);
    HTemp(EdgesToEstimate) = patt; % score this graph
    HTemp = HTemp + HTemp';
    ScoreBIC(k+1) = BICScore(X(:,cc),HTemp,HestSmall,gamma);
end

indMin = argmin(ScoreBIC);
patt = BinaryVector(indMin-1,NumEst);
HTemp = spalloc(length(cc),length(cc),length(cc));
HTemp(EdgesToEstimate) = patt;
HTemp = HTemp + HTemp';
p = size(X,2);
G = spalloc(p,p,p);
G(cc,cc) = HTemp;

end

function b = BinaryVector(d,N)

s = dec2bin(d,N);
b = zeros(1,N);
for k = 1:N
    b(k) = str2double(s(k));
end

end