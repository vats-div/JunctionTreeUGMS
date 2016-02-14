function [bhat, RetOptions] = myLasso(y,X,options)

% options.ReturnIndex = 1 -> bhat is an index set, else bhat is a vector

p = size(X,2);
n = size(X,1);

if p < n
    lse = inv(X'*X) * X' * y;
else
    lse = ones(p,1);
end

ReturnIndex = 1;
lambdaPattern = ones(p,1);
NumRegul = options.NumRegul;
SearchFurther = 1;

if isfield(options,'ReturnIndex')
    ReturnIndex = options.ReturnIndex;
end

if isfield(options,'lambdaPattern')
    lambdaPattern = options.lambdaPattern;
end

if isfield(options,'gamma')
    gamma = options.gamma;
end

if isfield(options,'SearchFurther')
    SearchFurther = options.SearchFurther;
end

%lambdaPattern = lambdaPattern .* 1./abs(lse).^2;

NumEst = double(sum(lambdaPattern>0));

if NumEst <= 5
%     % list out all possible options for lambdaPattern and compute BICScore
%     % for each option
%     gamma = 0;
%     indEst = (lambdaPattern > 0); % location of entries to estimate
%     BICScore = zeros(2^NumEst,1);
%     for k = 0:2^NumEst-1
%         patt = BinaryVector(k,NumEst);
%         supp = ones(p,1);
%         supp(indEst) = patt;
%         bhat{k+1} = EstimateBhat(y,X,supp);
%         sse = norm(y - X * bhat{k+1})^2;
%         NumNonZeros = sum(patt(:));
%         BICScore(k+1) = n*log(sse) + NumNonZeros*log(n) + 2*NumNonZeros*gamma*log(p);
%     end
%     indMin = argmin(BICScore);
%     bhatTemp = bhat{indMin}; clear bhat;
%     bhat = bhatTemp; clear bhatTemp;

XX = [y X];
Hest = zeros(p+1,p+1);
Hest(1,:) = 1; Hest(:,1) = 1; Hest = setdiag(Hest,0);
H = ones(p+1,p+1);
cc = 1:p+1;
Ghat = full(CheckEdgesRobust(XX,H,Hest,cc));
bhat = Ghat(1,2:end);
bhat = reshape(bhat,p,1);

else
    if options.METHOD == 1
        bhat = myLassoMainFunction(y,X,lambdaPattern,gamma,NumRegul,SearchFurther);
    else
        bhat = myLassoMainFunctionCV(y,X,lambdaPattern,gamma,NumRegul,SearchFurther);
    end
end

RetOptions.bhat = bhat;

if ReturnIndex == 1
    bhatPattern = bhat .* lambdaPattern;
    bhat = find(abs(bhatPattern) > 0);
end

end

function bhat = EstimateBhat(y,X,supp)

if (sum(supp(:)) == 0)
    bhat = supp;
    return;
else
    ne = find(supp ~= 0);
    bhat = zeros(size(X,2),1);
    bb = inv(X(:,ne)'*X(:,ne))*X(:,ne)' * y;
    bhat(ne) = bb;
end

end

function b = BinaryVector(d,N)

s = dec2bin(d,N);
b = zeros(1,N);
for k = 1:N
    b(k) = str2double(s(k));
end

end

function bhat = myLassoMainFunction(y,X,lambdaPattern,gamma,NumRegul,SearchFurther)

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

end

function [bhat, lambda] = myLassoMainFunctionCV(y,X,lambdaPattern,gamma,NumRegul,SearchFurther)

op = struct('verbose',0);
p = size(X,2);
n = size(X,1);
w = zeros(p,1);
[f,g] = SquaredError(w,X,y);
lambdaMax = max(abs(g));
lambdaMin = FindLambdaMin(y,X);
lambdaValues = linspace(lambdaMin,lambdaMax,NumRegul);
% if isfield(options,'N');
%     N = options.N;
% else
%     N = 30;
% end
N = 30;

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
	%display(['Fold = ' num2str(nf)]);
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
%funObj = @(b)SquaredError(b,X,y);
funObj = @(b)SquaredError(b,XTrain,yTrain);
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
