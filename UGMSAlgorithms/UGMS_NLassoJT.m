function [G, ReturnOptions] = UGMS_NLassoJT(X,H,options)

% Find junction tree fn H to find a graph
% ReturnOptions

UGMSAlg = @(X,H,Hest,options) UGMS_NLasso(X,H,Hest,options);
[G, ReturnOptions] = UGMS_RG(H,UGMSAlg,X,options);

end
