function G = UGMS_GLassoJT(X,H,options)

% junction tree based graphical Lasso
% H is the graph that contains (most of) the true edges.

UGMSAlg = @(X,H,Hest,options) UGMS_GLasso(X,H,Hest,options);
G = UGMS_RG(H,UGMSAlg,X,options);
G = sparse(G);


end
