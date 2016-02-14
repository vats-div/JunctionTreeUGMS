function SupportEstimate = ScreeningGroupLassoGivenLambda(y,X,clus_size,lambda,groups)

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

wt = reshape(wt,length(wt),1);
wt = ones(length(wt),1);

funObj = @(b)SquaredError(b,X,y);
binit = zeros(p,1);
bhat = L1GeneralGroup_Auxiliary(funObj,binit,lambda*wt,groups,options);
SupportEstimate = find(abs(bhat) > 0);

end
