function [invSigma, G] = NeighborhoodGraph(p,d,ss)

% function G = NeighborhoodGraph(p,d,ss)
% modifies from geo in contest
% generate points on [0,1]^2
% G(i,j) = 1/d-e with probability (2sqrt(pi))^-1 * exp(-4 * norm(i-j)^2)
% else zero
% if ss exists, then each entry in the inverse covariance is random
% restrict the degree of the graph to be less than d

m = 2;
coords = rand(p,m);
I = [];
J = [];
rho = 1/d - 0.005;

if nargin == 3
    RhoMat = rho * rand(p,p) .* randomSign(p,p);
    RhoMat = triu(RhoMat) + triu(RhoMat)';
else
    RhoMat = rho * ones(p,p);
end

for i = 2:p
    for j = 1:(i-1)
        diff = abs(coords(i,:) - coords(j,:));
        pr = 1/sqrt(2*pi) * exp(-4*norm(diff)^2);
        bool = rand(1) < pr;
        if bool == 1
            J = cat(1,J,j);
            I = cat(1,I,i);
        end
    end
end

S = ones(length(I),1);
G = sparse([I;J],[J;I],[S;S],p,p);
G = setdiag(G,0);

deg = sum(G);
inddeg = find(deg > d);

%Gtrue = zeros(p,p);

for k = 1:length(inddeg)
    node = inddeg(k);
    neigh = FindNeighborsUndirected(G,node,1);
    rndnum = randperm(d+1);
    rndnum = d+1;
    neighextra = neigh(rndnum(1):end);
    G(node,neighextra) = 0;
    G(neighextra,node) = 0;
end

invSigma = G .* RhoMat;
invSigma = setdiag(invSigma,1);

end
