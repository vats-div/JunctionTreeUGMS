function G = UGMS_PCJT(X,H,options)

% options.xi -> threshold
% options.jtType -> Type of junction tree
% options.pcType -> 1 for marginal and 2 for neighbors
% options.G -> graph to start from

UGMSAlg = @(X,H,Hest,options) UGMS_PC(X,H,Hest,options);
G = UGMS_RG(H,UGMSAlg,X,options);

end