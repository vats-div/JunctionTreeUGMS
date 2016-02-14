function SupportEstimate = ScreeningGroupLassoFast(y,X,clus_size,kE,groups)

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

% start with lambda large and then increment
bhat = zeros(p,1);
[f,g] = SquaredError(bhat,X,y);
lambda = max(abs(g));
q = kE;
bhat = L1GeneralGroup_Auxiliary(funObj,zeros(p,1),lambda*wt,groups,options);
numsupp = FindActiveNodesGroups(bhat,groups);

while numsupp < q
    lambda = lambda/2;
    bhat = L1GeneralGroup_Auxiliary(funObj,zeros(p,1),lambda,groups,options);
    numsupp = FindActiveNodesGroups(bhat,groups);
end

SupportEstimate = find(abs(bhat) > 0);

end


function numsupp = FindActiveNodesGroups(x_SpaRSA,groups)

% find the number of active groups in x_SpaRSA

x_SpaRSA = double(x_SpaRSA ~= 0);
x_SpaRSA = unique(x_SpaRSA .* groups);
x_SpaRSA = setdiff(x_SpaRSA,0);
numsupp = length(x_SpaRSA);

end
