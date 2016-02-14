function SupportEstimate = ScreeningGroupLasso(y,X,clus_size,kE,groups)

n = size(X,1);
p = size(X,2);
options.norm = 2;
options.verbose = 0;

% assign groups randomly between 1 and clus_size

if nargin == 4
display('finding groups')
num_nodes_assigned = 0;
indgroup = 0;
gg = [];
while num_nodes_assigned < p
    indgroup = indgroup + 1;
    num_nodes_group = randperm(clus_size);
    num_nodes_group = num_nodes_group(1);
    if num_nodes_assigned + num_nodes_group > p
        num_nodes_group = p - num_nodes_assigned;
    end
    wt(indgroup) = sqrt(num_nodes_group);
    gg = [gg indgroup * ones(1,num_nodes_group)];
    num_nodes_assigned = length(gg);
end

indd = randperm(p);
groups = gg(indd)';
else
	for k = 1:max(groups);
		wt(k) = sqrt(sum(groups == k));
	end
end

num_group = max(groups);
wt = reshape(wt,length(wt),1);

if (kE > num_group)
    SupportEstimate = 1:p;
    return;
end

funObj = @(b)SquaredError(b,X,y);

% find the lambda for which estimates are all zero first
% this is lambdaMax
binit = zeros(p,1);
numsupp = p;
q = kE;
tau = 0.01;

while numsupp > q
        lambda = tau * (sqrt(log(p)/n)) * n * 0.5;
        bhat = L1GeneralGroup_Auxiliary(funObj,binit,lambda*wt,groups,options);
        numsupp = FindActiveNodesGroups(bhat,groups);
        tau = tau*2;
end

lambdaMax = tau/1.2;

% Now find lambaMin

while numsupp < q
        lambda = tau * (sqrt(log(p)/n)) * n * 0.5;
        bhat = L1GeneralGroup_Auxiliary(funObj,binit,lambda*wt,groups,options);
        numsupp = FindActiveNodesGroups(bhat,groups);
        tau = tau/2;
end

lambdaMin = tau*1.2;

lambdaValues = linspace(lambdaMax,lambdaMin,25) * (sqrt(log(p)/n)) * n * 0.5;

numAbove = 0;
indlambda = 0;
for lambda = lambdaValues
	indlambda = indlambda + 1;
        bhat = L1GeneralGroup_Auxiliary(funObj,binit,lambda*wt,groups,options);
        numsupp(indlambda) = FindActiveNodesGroups(bhat,groups);
	if numsupp(indlambda) == q
		SupportEstimate = find(abs(bhat) > 0);
		return;
	end
	if numsupp(indlambda) > q
		numAbove = numAbove + 1;
	end
	if numAbove > 3
		break;
	end
end

% find the lambda for which numsupp is smallest above q

numsupp(numsupp < q) = q * 1000;
[temp, indd] = min(numsupp);
lambda = lambdaValues(indd);
bhat = L1GeneralGroup_Auxiliary(funObj,binit,lambda*wt,groups,options);
SupportEstimate = find(abs(bhat) > 0);

end


function numsupp = FindActiveNodesGroups(x_SpaRSA,groups)

% find the number of active groups in x_SpaRSA

x_SpaRSA = double(x_SpaRSA ~= 0);
x_SpaRSA = unique(x_SpaRSA .* groups);
x_SpaRSA = setdiff(x_SpaRSA,0);
numsupp = length(x_SpaRSA);

end
