% example to check whether all functions are working

% generate a Gaussian graphical model

p = 100; % number of vertices
n = 100;
kappa = 1;
perc = 0.2; % percentage of weak edges
p1 = round(perc * p);
d1 = 8; rho1 = 0.15;
d2 = 4; rho2 = 1/4-0.001;
[invSigma, Gtrue] = TwoNeighborhoodGraphs(p,d1,d2,perc,rho1,rho2);
Sigma = inv(invSigma);

% Gtrue is the graph
% Sigma is covariance matrix

% generate observations
X = mvnrnd(zeros(p,1),Sigma,n);


% options for screening algorithm
opScr.dmax = 0; opScr.K = 0; opScr.m = 1; % not important for the script
opScr.eta = kappa; opScr.H = ones(p,p);
opScr.numSearch = 5; opScr.numRepeat = 2;

% set options for using GLasso (graphical Lasso)
opGLasso.EstimateLambda = 1;
opGLasso.SearchLambda = linspace(0.01,1.0,10);
opGLasso.Nodes = 1:p;
opGLasso.Gtrue = Gtrue;
opGLasso.gamma = 0.5;
opGLasso.SearchFurther = 1;
opGLasso.jtType = 2;
opGLasso.displayInd = 0;
opGLasso.SepSize = kappa+1;

% set options for using NLasso (neighborhood selection using Lasso)
opNLasso.NumRegul = 10;
opNLasso.gamma = 0.5;
opNLasso.SearchFurther = 1;
opNLasso.Nodes = 1:p;
opNLasso.jtType = 2;
opNLasso.displayInd = 0;
opNLasso.SepSize = kappa+1;

% set options for UGMS_PC (PC-Algorithm)
opPC.pcType = 1;
opPC.EstimateXi = 1;
opPC.SearchFurther = 1;
opPC.jtType = 2;
opPC.Delta = 0:kappa+1;
opPC.Gtrue = Gtrue;
opPC.displayInd = 0;
opPC.SepSize = kappa+1;


% First perform screen to compute H

display('Starting Graph Screening');
H = GraphScreening(X,opScr); H = setdiag(H,0);
display('Done Graph Screening');


% compute graphical Lasso estimates
% without using junction trees
opGLasso.METHOD = 4;
opGLasso.gamma = 0.5;
display('Running graphical Lasso');
Ghat_NOJT_EBIC_GL = UGMS_GLasso(X,H,H,opGLasso);
display('Done running graphical Lasso');

% JT+EBIC+GLasso
opGLasso.METHOD = 4;
display('Running junction tree based graphical Lasso');
Ghat_JT_EBIC_GL = UGMS_GLassoJT(X,H,opGLasso);
display('Done running junction tree based graphical lasso');

% PC-Algorithm without using junction trees
opPC.SearchXi = linspace(0.05,0.4,10);
opPC.METHOD = 4; % for EBIC
opPC.gamma = 0.5;
display('Running PC-Algorithm');
Ghat_NOJT_EBIC_PC = UGMS_PC(X,H,H,opPC);
display('Done running PC-Algorithm');

% JT+EBIC+PC
opPC.METHOD = 4;
opPC.gamma = 0.5;
opPC.SearchXi = linspace(0.05,0.4,10);
display('Running junction tree based PC-Algorithm');
Ghat_JT_EBIC_PC = UGMS_PCJT(X,H,opPC);
display('Done running junction tree based PC-Algorithm');

display('Evaluating Algorithms')

% WEDR is weak edge discovery rate
display(' ');
display(['WEDR - gLasso: ' num2str(TruePositiveGraphPartial(Gtrue,Ghat_NOJT_EBIC_GL,1,p1))]);
display(['TPR - gLasso: ' num2str(TruePositiveGraph(Gtrue,Ghat_NOJT_EBIC_GL,1))]);
display(['FPR - gLasso: ' num2str(FalsePositiveGraph(Gtrue,Ghat_NOJT_EBIC_GL,1))]);
display(['ED - gLasso: ' num2str(CompareGraphs(Gtrue,Ghat_NOJT_EBIC_GL))]);
display(['Number of edges -  gLasso: ' num2str(sum(Ghat_NOJT_EBIC_GL(:))/2)]);

display(' ');

display(['WEDR - JT-gLasso: ' num2str(TruePositiveGraphPartial(Gtrue,Ghat_JT_EBIC_GL,1,p1))]);
display(['TPR - JT-gLasso: ' num2str(TruePositiveGraph(Gtrue,Ghat_JT_EBIC_GL,1))]);
display(['FPR - JT-gLasso: ' num2str(FalsePositiveGraph(Gtrue,Ghat_JT_EBIC_GL,1))]);
display(['ED - JT-gLasso: ' num2str(CompareGraphs(Gtrue,Ghat_JT_EBIC_GL))]);
display(['Number of edges -  JT-gLasso: ' num2str(sum(Ghat_JT_EBIC_GL(:))/2)]);

% WEDR is weak edge discovery rate
display(' ');
display(['WEDR - PC: ' num2str(TruePositiveGraphPartial(Gtrue,Ghat_NOJT_EBIC_PC,1,p1))]);
display(['TPR - PC: ' num2str(TruePositiveGraph(Gtrue,Ghat_NOJT_EBIC_PC,1))]);
display(['FPR - PC: ' num2str(FalsePositiveGraph(Gtrue,Ghat_NOJT_EBIC_PC,1))]);
display(['ED - PC: ' num2str(CompareGraphs(Gtrue,Ghat_NOJT_EBIC_PC))]);
display(['Number of edges -  PC: ' num2str(sum(Ghat_NOJT_EBIC_PC(:))/2)]);

display(' ');

display(['WEDR - JT-PC: ' num2str(TruePositiveGraphPartial(Gtrue,Ghat_JT_EBIC_PC,1,p1))]);
display(['TPR - JT-PC: ' num2str(TruePositiveGraph(Gtrue,Ghat_JT_EBIC_PC,1))]);
display(['FPR - JT-PC: ' num2str(FalsePositiveGraph(Gtrue,Ghat_JT_EBIC_PC,1))]);
display(['ED - JT-PC: ' num2str(CompareGraphs(Gtrue,Ghat_JT_EBIC_PC))]);
display(['Number of edges -  JT-PC: ' num2str(sum(Ghat_JT_EBIC_PC(:))/2)]);