function [invSigma, Gtrue] = TwoChainGraphs(p,perc,rho1,rho2)

% function [invSigma, Gtrue] = TwoChainGraphs(p,ratio,rho1,rho2)
% Chain Graph
% p - number of nodes
% perc - percentage of smaller nodes
% rho1 - minimum partial correlation over smaller nodes
% rho2 - minimum partial correlation over larger nodes

p1 = perc*p; 
p2 = 1 - p1;

Gtrue = diag(ones(p-1,1),1) + diag(ones(p-1,1),-1);

invSigma = Gtrue * rho2;
invSigma(1:p1,1:p1) = Gtrue(1:p1,1:p1) * rho1;
invSigma = invSigma + eye(p,p);

end