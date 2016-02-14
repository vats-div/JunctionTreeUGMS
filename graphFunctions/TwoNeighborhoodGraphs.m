function [invSigma, Gtrue] = TwoNeighborhoodGraphs(p,d1,d2,perc,rho1,rho2)

% function [invSigma, Gtrue] = TwoNeighborhoodGraphs(p,perc,d1,d2,rho1,rho2)

% p: Number of nodes
% perc: percentage of smaller set
% d1: degree of smaller set (should be larger than d2)
% d2: degree of larger set
% rho1: partial correlation of smaller set
% rho2: partial correlation of larger set

p1 = ceil(perc * p);
if p1 < 4
    p1 = 4;
end
p2 = p - p1;
[invSigma1, Gtrue1] = NeighborhoodGraph(p1,d1); % can also use 
[invSigma2, Gtrue2] = NeighborhoodGraph(p2,d2); % can also use

% connect the last d1 nodes in Gtrue1 to the first d2 nodes in Gtrue2

Gtrue = zeros(p,p);
Gtrue(1:p1,1:p1) = Gtrue1;
Gtrue(p1+1:end,p1+1:end) = Gtrue2;

Gtrue(p1-2,p1+1) = 1; Gtrue(p1+1,p1-2) = 1;
Gtrue(p1-1,p1+2) = 1; Gtrue(p1+2,p1-1) = 1;
Gtrue(p1,p1+3) = 1; Gtrue(p1+3,p1) = 1;
Gtrue(p1-2,p1+2) = 1; Gtrue(p1+2,p1-2) = 1;

invSigma = zeros(p,p);
invSigma(1:p1,1:p1) = rho1 * Gtrue1;
invSigma(p1+1:end,p1+1:end) = rho2 * Gtrue2;
invSigma(p1-2,p1+1) = rho1; invSigma(p1+1,p1-2) = rho1;
invSigma(p1-1,p1+2) = rho1; invSigma(p1+2,p1-1) = rho1;
invSigma(p1,p1+3) = rho1; invSigma(p1+3,p1) = rho1;
invSigma(p1-2,p1+2) = rho1; invSigma(p1+2,p1-2) = rho1;

invSigma = invSigma + eye(p,p);

end