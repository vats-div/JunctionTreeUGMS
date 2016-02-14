function H = GraphScreening(X,options)

% options.dmax -> upper bound on maximum degree in the graph
% options.K -> Number of times to run Multiple grouping
% options.m -> Group size
% options.eta -> size of separators to search over
% options.numSearch -> Number of separators to search over
% options.numRepeat -> Number of times to repeat the search
% options.isGaussian -> 1 if Gaussian distribution
% options.H -> if H is given do not apply GhatMuG

% Read options
[dmax,K,m,eta,numSearch,numRepeat,isGaussian] = getAlloptions(options);

if isGaussian == 0
    display('GraphScreening only works for Gaussian distributions');
    H = [];
    return;
end

if ~isfield(options,'H')
    % Perform multiple grouping for screening
    H = MuGScreeningGraphsFast(X,K,m,dmax);
else
    H = options.H;
end

clear options;

for kk = 1:numRepeat
    for k = 0:eta
        H = CheckEdgesCrossValidation(X,H,k,numSearch);
    end
end

display('Done Screening');
    
end

function [dmax,K,m,eta,numSearch,numRepeat,isGaussian] = getAlloptions(options)

% deafault values
K = 10; m = 2; eta = 1; numSearch = 1; numRepeat = 1; isGaussian = 1;

if isfield(options,'dmax')
    dmax = options.dmax;
end

if isfield(options,'K')
    K = options.K;
end
if isfield(options,'m');
    m = options.m;
end
if isfield(options,'eta');
    eta = options.eta;
end
if isfield(options,'numSearch')
    numSearch = options.numSearch;
end
if isfield(options,'numRepeat')
    numRepeat = options.numRepeat;
end
if isfield(options,'isGaussian')
    isGaussian = options.isGaussian;
end

end