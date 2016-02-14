function [SupportEstimate,supp] = MuGScreeningGroupLassoGivenLambda(y,X,K,m,lambda,SupportEstimate)

n = size(X,1);
p = size(X,2);
options.norm = 2;
options.verbose = 0;

if nargin == 5
    SupportEstimate = 1:p;
    includesupp = 1;
else
    includesupp = 0;
end
   
groupChoice = 1; % outer 0 for inner

for indk = 1:K
	display(['MuG with K = ' num2str(indk)]);
	if groupChoice == 1
		groups = FindGroupsOuter(SupportEstimate,p,m);
	else
		groups = FindGroupsInner(SupportEstimate,p,m);
	end
	ss = ScreeningGroupLassoGivenLambda(y,X,m,lambda,groups);
    if indk == 1 && includesupp == 0
        SupposeEstimate = ss;
    else
        SupportEstimate = intersect(SupportEstimate,ss);
    end
	length(SupportEstimate)
	supp{indk+1} = SupportEstimate;
end


end

function groups = FindGroupsOuter(ss,p,m)

% Associate each element in ss to m-1 elements in 1:p \ ss

ss = reshape(ss,length(ss),1);
ssdiff = setdiff(1:p,ss);
ssdiff = reshape(ssdiff,length(ssdiff),1);
groups = zeros(p,1);
indgroup = 0;
while (isempty(ss) == 0)
	sz = randperm(m);
    sz = m;
	sz = sz(1);
	sz = min([sz,length(ssdiff)]);
	indgroup = indgroup + 1;
	ssdiff = ssdiff(randperm(length(ssdiff)));
	ggind = [ss(1)  ; ssdiff(1:sz-1)];
	groups(ggind) = indgroup;
	ss(1) = [];
	ssdiff(1:sz-1) = [];
	if isempty(ss)
		break;
	end
	if isempty(ssdiff)
		break;
	end
end

ssLeft = [ss ; ssdiff]; % nodes left

ssLeft = ssLeft(randperm(length(ssLeft)));
while (isempty(ssLeft) == 0)
	if length(ssLeft) < m
		groups(ssLeft) = indgroup + 1;
		break;
	end
	sz = randperm(m);
	sz = sz(1);
	indgroup = indgroup + 1;
	ggind = [ssLeft(1:sz)];
	ssLeft(1:sz) = [];
	groups(ggind) = indgroup;
	if length(ssLeft) < m
		groups(ssLeft) = indgroup + 1;
		break;
	end
end

end

function groups = FindGroupsInner(ss,p,m)

% group elements in ss together and then group elements in ssdiff together

ss = reshape(ss,length(ss),1);
ssdiff = setdiff(1:p,ss);
ssdiff = reshape(ssdiff,length(ssdiff),1);
groups = zeros(p,1);
indgroup = 0;

while (isempty(ss) == 0)
        indgroup = indgroup + 1;
        ssdiff = ssdiff(randperm(length(ssdiff)));
        ggind = ss(1:m);
        groups(ggind) = indgroup;
        ss(1:m) = [];
        if isempty(ss)
                break;
        end
        if length(ss) < m
                groups(ss) = indgroup + 1;
                break;
        end
	ssAll = [ss ; ssdiff];
	if (indgroup + length(ssAll) < n + 10)
		ssAll = [ss ; ssdiff];
		for ll = 1:length(ssAll)
			indgroup = indgroup + 1;
			groups(ssAll(ll)) = indgroup;
		end
		return;
	end
end

while (isempty(ssdiff) == 0)
        indgroup = indgroup + 1;
        ggind = [ssdiff(1:m)];
        ssdiff(1:m) = [];
        groups(ggind) = indgroup;
        if length(ssdiff) < m
                groups(ssdiff) = indgroup + 1;
                break;
        end
	ssAll = ssdiff;
	if (indgroup + length(ssAll) < n + 10)
		ssAll = [ss ; ssdiff];
		for ll = 1:length(ssAll)
			indgroup = indgroup + 1;
			groups(ssAll(ll)) = indgroup;
		end
		return;
	end
end

end
