% Parameters for specifying graph
p = 300; % number of nodes
perc = 0.1;
onem = ones(p,p);
display('NG+NG')
N = [300 500 700];

d1 = 8; rho1 = 0.15;
d2 = 4; rho2 = 1/4-0.001;
[invSigma, Gtrue] = TwoNeighborhoodGraphs(p,d1,d2,perc,rho1,rho2);
save GtrueNGNGp300 invSigma Gtrue;
Sigma = inv(invSigma);
savefile = 'SimulationsData/NGNGModelSelectionp300METHOD2';
ExpUGMSNGNGp300;
%clear all;

%load GtrueNGNGp300;
% Parameters for specifying graph
%p = 300; % number of nodes
%perc = 0.1;
%p1 = perc * p;
%onem = ones(p,p);
%display('NG+NG')
%N = [200 300 500 700 1000 ];
%invSigma(1:p1,1:p1) = 0.075;
%invSigma = setdiag(invSigma,1);
%Sigma = inv(invSigma);
%savefile = 'SimulationsData/NGNGModelSelectionp100METHOD3';
%ExpUGMSNGNGp300;
