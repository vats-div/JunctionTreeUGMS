% RunAllExperiments

%% NON JUNCTION TREE BASED METHOD

% Ghat_NOJT_Oracle_GL
% Ghat_NOJT_EBIC_GL
% Ghat_NOJT_EBIC_PC
% Ghat_NOJT_Oracle_PC
% Ghat_NOJT_EBIC25_NL
% Ghat_NOJT_EBIC50_NL
% Ghat_NOJT_EBIC75_NL
% Ghat_NOJT_EBIC100_NL
% Ghat_JT_EBIC_GL
% Ghat_NOJT_EDGES_GL
% Ghat_JT_EBIC25_NL
% Ghat_JT_EBIC50_NL
% Ghat_JT_EBIC75_NL
% Ghat_JT_EBIC100_NL
% Ghat_JT_EBIC_PC
% Ghat_NOJT_EDGES_PC

%% SCREENING ALGORITHM

opScr.H = ones(p,p); % this way, the Multiple grouping screening approach is not used
H = GraphScreening(X,opScr);

display('Done Graph Screening');

%% NON JUNCTION TREE BASED METHOD

% Oracle method that uses the true graph and GLasso to find an estimate
opGLasso.METHOD = 3;
Ghat_NOJT_Oracle_GL = UGMS_GLasso(X,H,H,opGLasso);
display('Computed Oracle Graph using GLasso');

% Extended BIC and GLasso
opGLasso.METHOD = 4;
opGLasso.gamma = 0.5;
Ghat_NOJT_EBIC_GL = UGMS_GLasso(X,H,H,opGLasso);
display('Computed EBIC Graph using GLasso');

% Extended BIC and PC-Algorithm
opPC.SearchXi = linspace(0.05,0.4,10);
opPC.METHOD = 4; % for EBIC
opPC.gamma = 0.5;
Ghat_NOJT_EBIC_PC = UGMS_PC(X,H,H,opPC);
display('Done with EBIC+PC');

opPC.METHOD = 3;
Ghat_NOJT_Oracle_PC = UGMS_PC(X,H,H,opPC);

% Extended BIC and NLasso
opNLasso.METHOD = 1;
opNLasso.gamma = 0.5;
Ghat_NOJT_EBIC25_NL = UGMS_NLasso(X,H,H,opNLasso);
opNLasso.gamma = 0.5;
Ghat_NOJT_EBIC50_NL = UGMS_NLasso(X,H,H,opNLasso);
opNLasso.gamma = 0.75;
Ghat_NOJT_EBIC75_NL = UGMS_NLasso(X,H,H,opNLasso);
opNLasso.gamma = 1.00;
Ghat_NOJT_EBIC100_NL = UGMS_NLasso(X,H,H,opNLasso);

display('Done with EBIC+NLasso');

%% JUNCTION TREE BASED APPROACHES

% JT+EBIC+GLasso
opGLasso.METHOD = 4;
Ghat_JT_EBIC_GL = UGMS_GLassoJT(X,H,opGLasso);
display('Done with JT+EBIC+GLasso');

% Compute GL graph with same number of edges
%opGLasso.TP = TruePositiveGraph(Gtrue,Ghat_JT_EBIC_GL,1);
opGLasso.METHOD = 5;
opGLasso.NumEdges = sum(Ghat_JT_EBIC_GL(:))/2;
Ghat_NOJT_EDGES_GL = UGMS_GLasso(X,H,H,opGLasso);

% JT+EBIC+NLasso
opNLasso.METHOD = 1;
opNLasso.gamma = 0.25;
Ghat_JT_EBIC25_NL = UGMS_NLassoJT(X,H,opNLasso);
opNLasso.gamma = 0.5;
Ghat_JT_EBIC50_NL = UGMS_NLassoJT(X,H,opNLasso);
opNLasso.gamma = 0.75;
Ghat_JT_EBIC75_NL = UGMS_NLassoJT(X,H,opNLasso);
opNLasso.gamma = 1.00;
Ghat_JT_EBIC100_NL = UGMS_NLassoJT(X,H,opNLasso);
display('Done with JT+EBIC+NLasso');

% JT+EBIC+PC
opPC.METHOD = 4;
opPC.gamma = 0.5;
opPC.SearchXi = linspace(0.05,0.4,10);
Ghat_JT_EBIC_PC = UGMS_PCJT(X,H,opPC);
display('Done with JT+EBIC+PC');

% Compute PC graph with same EDGES
opPC.METHOD = 5;
opPC.NumEdges = sum(Ghat_JT_EBIC_PC(:))/2;
Ghat_NOJT_EDGES_PC = UGMS_PC(X,H,H,opPC);
