function [invSigma, Gtrue] = TwoCycleGraphs(p,perc,rho1,rho2)

% function [invSigma, Gtrue] = TwoChainGraphs(p,ratio,rho1,rho2)
% Chain Graph
% p - number of nodes
% perc - percentage of smaller nodes
% rho1 - minimum partial correlation over smaller nodes
% rho2 - minimum partial correlation over larger nodes

p1 = ceil(perc*p); 

% cylcle of length 4
Gtrue = diag(ones(p-1,1),1) + diag(ones(p-1,1),-1);
for k = 1:2:p1
    Gtrue(k,k+2) = 1;
    Gtrue(k+2,k) = 1;
end
for k = 2:2:p1
    Gtrue(k,k+2) = 1;
    Gtrue(k+2,k) = 1;
end


Gtrue(1,p1) = 1; Gtrue(p1,1) = 1;
Gtrue(p1+1,end) = 1; Gtrue(end,p1+1) = 1;

for k = p1:5:p-1
    Gtrue(k,min(p,k+5)) = 1;
    Gtrue(min(p,k+5),k) = 1;
end

invSigma = Gtrue * rho2;
invSigma(1:p1,1:p1) = Gtrue(1:p1,1:p1) * rho1;
invSigma = invSigma + eye(p,p);

end