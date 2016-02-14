function [bhatRet, returnOptions] = LassoWithModelSelection(y,X,options)

% function [bhat, returnOptions] = LassoWithModelSelection(y,X,options)

% options.lambda -> if lambda is a scalar, then regularization is \lambda
% ||b||_1, else a weighted l1 norm

% options.METHOD -> 0 for using given lambda, 1 for BIC, 2 for stability
% selection, 3 for support size, 4 for CV
% options.btrue -> optimal value is known
% options.lambda -> given lambda

% change X to include intercept
%Xc = zeros(size(X,1),size(X,2)+1);
%Xc(:,2:end) = X;
%Xc(1,1) = 1;
%X = Xc; clear Xc;

if nargout == 2
    returnOptions = [];
end

if nargin == 2
    options.METHOD = 0;
end

% if ~isa(fsX,'function_handle')
%     % make it a function handle
%     X = fsX; clear fsX;
%     fsX = @(w,mode) fun_X(w,mode,X,X');
% end

pb = size(X,2);
ReturnIndex = 1;
lambdaPattern = ones(pb,1);
SearchFurther = 1;
gamma = 0;
useMethod = 0;
lambda = 1.0;
NumRegul = 50;
indNonZeros = (ones(pb,1) > 0);

if isfield(options,'lambda')
    lambda = options.lambda;
end

if isfield(options,'k')
    SuppSize = options.k;
end

if isfield(options,'ReturnIndex')
    ReturnIndex = options.ReturnIndex;
end

if isfield(options,'lambdaPattern')
    lambdaPattern = options.lambdaPattern;
    indNonZeros = (~isinf(lambdaPattern));
    X = X(:,indNonZeros);
    lambdaPattern = lambdaPattern(indNonZeros);
end
if isfield(options,'gamma')
    gamma = options.gamma;
end
if isfield(options,'METHOD')
    useMethod = options.METHOD;
end

p = size(X,2);

switch useMethod
    case 0
        w = zeros(p,1);
        funObj = @(w)SquaredError(w,X,y);
        op = struct('verbose',0);
        bhat = L1General2_PSSgb(funObj,w,lambda*lambdaPattern,op);
    case 1
        [bhat, lambda] = myLassoMainFunction(y,X,lambdaPattern,gamma,NumRegul,SearchFurther);
        returnOptions.lambda = lambda;
    case 2
        bhat = myLassoMainFunctionSTAB(y,X,lambdaPattern,NumRegul,options);
    case 3
        bhat = myLassoMainFunctionSuppSize(y,X,lambdaPattern,NumRegul,SuppSize);
    case 4
        [bhat, lambda] = myLassoMainFunctionCV(y,X,lambdaPattern,NumRegul,options);
        returnOptions.lambda = lambda;
    otherwise
        bhat = myLassoMainFunction(y,X,lambdaPattern,gamma,NumRegul,SearchFurther);
        
end

bhatRet = zeros(pb,1);
bhatRet(indNonZeros) = bhat;
end

function [bhat, lambda] = myLassoMainFunctionSuppSize(y,X,lambdaPattern,NumRegul,SuppSize)

op = struct('verbose',0);
p = size(X,2);
n = size(X,1);
w = zeros(p,1);
[f,g] = SquaredError(w,X,y);
lambdaMax = max(abs(g));
lambdaMin = 0.2/n;
C = (log(lambdaMax) - log(lambdaMin))/NumRegul;
lambdaValues = exp(-(1:NumRegul) * C) * lambdaMax;
lambdaValues = sort(lambdaValues,'ascend');

err = zeros(length(lambdaValues),1);

binit = zeros(p,1);
funObj = @(b)SquaredError(b,X,y);
indlambda = 0;
for lambda = lambdaValues
    indlambda = indlambda + 1;
    bhat(:,indlambda) = L1General2_PSSgb(funObj,binit,lambda*lambdaPattern,op);
    binit = bhat(:,indlambda);
    err(indlambda) = length(find(binit)) - SuppSize;
end

err(err < 0) = 100000;
[minerr indmin] = min(err);

funObj = @(b)SquaredError(b,X,y);
binit = zeros(p,1);
bhat = L1General2_PSSgb(funObj,binit,lambdaValues(indmin)*lambdaPattern,op);

ss = find(abs(bhat) > 0);
bhat = zeros(p,1);
if ~isempty(ss)
    Xss = X(:,ss);
    bhat(ss) = inv(Xss'*Xss) * Xss' * y;
end

end

function [bhat, lambda] = myLassoMainFunctionCV(y,X,lambdaPattern,NumRegul,options)

op = struct('verbose',0);
p = size(X,2);
n = size(X,1);
w = zeros(p,1);
[f,g] = SquaredError(w,X,y);
lambdaMax = max(abs(g));
lambdaMin = FindLambdaMin(y,X);
%lambdaValues = linspace(lambdaMin,lambdaMax,NumRegul);

C = (log(lambdaMax) - log(lambdaMin))/NumRegul;
lambdaValues = exp(-(1:NumRegul) * C) * lambdaMax;
lambdaValues = sort(lambdaValues,'ascend');

if isfield(options,'N');
    N = options.N;
else
    N = 30;
end

numTrain = round(0.7 * n);
for nf = 1:N
	cc = randperm(n);
	ccT{nf} = cc(1:numTrain);
	ccTe{nf} = cc(numTrain+1:end);
end
err = zeros(length(lambdaValues),1);
for nf = 1:N
	indTrain = ccT{nf};
 	indTest = ccTe{nf};
	XTrain = X(indTrain,:);
	yTrain = y(indTrain,:);
	XTest = X(indTest,:);
	yTest = y(indTest,:);
	display(['Fold = ' num2str(nf)]);
	% for each lambda, compute bhat on Training
	indlambda = 0;
	binit = zeros(p,1);
	funObj = @(b)SquaredError(b,XTrain,yTrain);
	for lambda = lambdaValues
		indlambda = indlambda + 1;
    		bhat(:,indlambda) = L1General2_PSSgb(funObj,binit,lambda*lambdaPattern,op);
		binit = bhat(:,indlambda);
		err(indlambda) = err(indlambda) + norm(yTest - XTest * bhat(:,indlambda))^2;
	end
end

[minerr indmin] = min(err);

% indmin is the index chosen
funObj = @(b)SquaredError(b,X,y);
%funObj = @(b)SquaredError(b,XTrain,yTrain);
binit = zeros(p,1);
bhat = L1General2_PSSgb(funObj,binit,lambdaValues(indmin)*lambdaPattern,op);

ss = find(abs(bhat) > 0);
bhat = zeros(p,1);
if ~isempty(ss)
    Xss = X(:,ss);
    bhat(ss) = inv(Xss'*Xss) * Xss' * y;
end

lambda = lambdaValues(indmin);

end

function bhat = myLassoMainFunctionSTAB(y,X,lambdaPattern,NumRegul,options)

p = size(X,2);
n = size(X,1);
w = zeros(p,1);
[f,g] = SquaredError(w,X,y);
lambdaMax = max(abs(g));
lambdaMin = FindLambdaMin(y,X);
SearchLambda = linspace(lambdaMin,lambdaMax,NumRegul);
thresh = options.thresh;
N = options.N;

PiSel = zeros(p,NumRegul);
bhat = zeros(p,1);
funObj = @(b)SquaredError(b,X,y);
binit = zeros(p,1);
op = struct('verbose',0);
b = floor(3*n/4);
q = n;
for indn = 1:N
	indT = randperm(n);
	indT = indT(1:b); 
	XTrain = X(indT,:);
	yTrain = y(indT,:);
	funObj = @(b)SquaredError(b,XTrain,yTrain);
 	indlambda = 0;
	binit = zeros(p,1);
	for lambda = SearchLambda 
		indlambda = indlambda + 1;
	        lambdaVec = lambda * lambdaPattern;
	        bb = L1General2_PSSgb(funObj,binit,lambdaVec,op);
            NumNon = length(find(abs(bb) > 0));
		if NumNon > q
			break;
		else
			PiSel(:,indlambda) = PiSel(:,indlambda) + (abs(bb) > 0);
		end
	end
end

PiSel = PiSel/N;
PiSel = max(PiSel,[],2);

% do stability selection

ss = find(PiSel > thresh);
bhat = zeros(p,1);
if ~isempty(ss)
    Xss = X(:,ss);
    bhat(ss) = inv(Xss'*Xss) * Xss' * y;
end

end

function [bhat, lambda] = myLassoMainFunction(y,X,lambdaPattern,gamma,NumRegul,SearchFurther)

p = size(X,2);
n = size(X,1);
w = zeros(p,1);
[f,g] = SquaredError(w,X,y);
lambdaMax = max(abs(g));
lambdaMin = FindLambdaMin(y,X);
SearchLambda = linspace(lambdaMin,lambdaMax,NumRegul);

indk = 0;
ScoreBIC = zeros(length(SearchLambda),1);

op = struct('verbose',0);

for lambda = SearchLambda
    funObj = @(w)SquaredError(w,X,y);
    bhat = L1General2_PSSgb(funObj,w,lambda.*lambdaPattern,op);
    indk = indk + 1;
    sse = norm(y - X * bhat)^2;
    NumNonZeros = sum(bhat .* lambdaPattern ~= 0);
    ScoreBIC(indk) = n*log(sse) + NumNonZeros*log(n) + 2*NumNonZeros*gamma*log(p);
end

indMin = argmin(ScoreBIC);

%NumRegul = round(NumRegul/2);

if SearchFurther == 1
    if indMin > 1 && indMin < length(SearchLambda)
        SearchLambda = linspace(SearchLambda(indMin-1),SearchLambda(indMin+1),NumRegul);
    else
        SearchLambda = linspace(SearchLambda(indMin)/2,SearchLambda(indMin)*2,NumRegul);
        
    end
    clear ScoreBIC;
    indk = 0;
    for lambda = SearchLambda
        funObj = @(w)SquaredError(w,X,y);
        bhat = L1General2_PSSgb(funObj,w,lambda.*lambdaPattern,op);
        indk = indk + 1;
        sse = norm(y - X * bhat)^2;
        NumNonZeros = sum(bhat .* lambdaPattern ~= 0);
        ScoreBIC(indk) = n*log(sse) + NumNonZeros*log(n) + 2*NumNonZeros*gamma*log(p);
    end
    indMin = argmin(ScoreBIC);
end

bhat = L1General2_PSSgb(funObj,w,SearchLambda(indMin)*lambdaPattern,op);
lambda = SearchLambda(indMin);

end

function lambda = FindLambdaMin(y,X)

p = size(X,2);
n = size(X,1);
minn = min(n,p);
options = struct('verbose',0);
w = zeros(p,1);
lambda = 1;
funObj = @(w)SquaredError(w,X,y);
bhat = L1General2_PSSgb(funObj,w,lambda*ones(p,1),options);
NumNonZeros = sum(bhat ~= 0); 
numTimes = 0;
while (NumNonZeros < minn && numTimes < 20)
    lambda = lambda/2;
    bhat = L1General2_PSSgb(funObj,w,lambda*ones(p,1),options);
    NumNonZeros = sum(bhat ~= 0); 
    numTimes = numTimes + 1;
end

end
