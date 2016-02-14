% Parameters for specifying graph
p = 100; % number of nodes
perc = 0.2;
onem = ones(p,p);
% Generate Graph
% fix the graph

% CHAIN
display('Chain')

N = [300 500];
[invSigma, Gtrue] = TwoChainGraphs(p,perc,rho1,rho2);
Sigma = inv(invSigma);

% Parameters for the experiment
numTrials = 50; % number of times to do the experiment

% SET DIFFERENT OPTIONS

% options for screening algorithms
opScr.dmax = 30; opScr.K = 10; opScr.m = 2; 
opScr.eta = 0; 
opScr.numSearch = 1; opScr.numRepeat = 5;

% set options for using GLasso (graphical Lasso)
opGLasso.EstimateLambda = 1;
opGLasso.SearchLambda = linspace(0.01,1.0,10);
opGLasso.Nodes = 1:p;
opGLasso.Gtrue = Gtrue;
opGLasso.gamma = 0.5;
opGLasso.SearchFurther = 1;
opGLasso.jtType = 2;
opGLasso.displayInd = 1;
opGLasso.SepSize = 2;

% set options for using NLasso (neighborhood selection using Lasso)
opNLasso.NumRegul = 10;
opNLasso.gamma = 0.5;
opNLasso.SearchFurther = 1;
opNLasso.Nodes = 1:p;
opNLasso.jtType = 2;
opNLasso.displayInd = 1;
opNLasso.SepSize = 2;

% set options for UGMS_PC (PC-Algorithm)
opPC.pcType = 1;
opPC.EstimateXi = 1;
opPC.SearchFurther = 1;
opPC.jtType = 2;
opPC.Delta = 0:1;
opPC.Gtrue = Gtrue;
opPC.displayInd = 1;
opPC.SepSize = 2;

for i = 1:numTrials
    indn = 0;
    for n = N
        indn = indn + 1;
        XX{i,indn} = mvnrnd(zeros(p,1),Sigma,n);
        X = XX{i,indn};
        % Run All METHODS
        RunAllExperiments;
        Save_H{i,indn} = sparse(H);
        Save_Ghat_NOJT_Oracle_GL{i,indn} = sparse(Ghat_NOJT_Oracle_GL);
        Save_Ghat_NOJT_EBIC_GL{i,indn} = sparse(Ghat_NOJT_EBIC_GL);
        Save_Ghat_NOJT_EBIC_PC{i,indn} = sparse(Ghat_NOJT_EBIC_PC);
        Save_Ghat_NOJT_Oracle_PC{i,indn} = sparse(Ghat_NOJT_Oracle_PC);
        Save_Ghat_NOJT_EBIC25_NL{i,indn} = sparse(Ghat_NOJT_EBIC25_NL);
        Save_Ghat_NOJT_EBIC50_NL{i,indn} = sparse(Ghat_NOJT_EBIC50_NL);
        Save_Ghat_NOJT_EBIC75_NL{i,indn} = sparse(Ghat_NOJT_EBIC75_NL);
        Save_Ghat_NOJT_EBIC100_NL{i,indn} = sparse(Ghat_NOJT_EBIC100_NL);
        Save_Ghat_JT_EBIC_GL{i,indn} = sparse(Ghat_JT_EBIC_GL);
        Save_Ghat_NOJT_EDGES_GL{i,indn} = sparse(Ghat_NOJT_EDGES_GL);
        Save_Ghat_JT_EBIC25_NL{i,indn} = sparse(Ghat_JT_EBIC25_NL);
        Save_Ghat_JT_EBIC50_NL{i,indn} = sparse(Ghat_JT_EBIC50_NL);
        Save_Ghat_JT_EBIC75_NL{i,indn} = sparse(Ghat_JT_EBIC75_NL);
        Save_Ghat_JT_EBIC100_NL{i,indn} = sparse(Ghat_JT_EBIC100_NL);
        Save_Ghat_JT_EBIC_PC{i,indn} = sparse(Ghat_JT_EBIC_PC);
        Save_Ghat_NOJT_EDGES_PC{i,indn} = sparse(Ghat_NOJT_EDGES_PC);
        save(savefile);
    end
end
