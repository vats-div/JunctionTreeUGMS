function [SupportEstimate,supp] = MuGScreeningGroupLassoFast(y,X,K,m,kE)

n = size(X,1);
p = size(X,2);
options.norm = 2;
options.verbose = 0;

groups = 1:p';

if nargin == 4
	kE = n;
end

%ss = ScreeningLasso(y,X,kE);
ss = ScreeningGroupLassoFast(y,X,1,kE,(1:p)');
SupportEstimate = ss;
supp{1} = ss;

groupChoice = 1; % outer 0 for inner

for indk = 1:K
	display(['MuG with K = ' num2str(indk)]);
	if groupChoice == 1
		groups = FindGroupsOuter(SupportEstimate,p,m);
	else
		groups = FindGroupsInner(SupportEstimate,p,m);
	end
	ss = ScreeningGroupLassoFast(y,X,m,kE,groups);
	SupportEstimate = intersect(SupportEstimate,ss);
	supp{indk+1} = SupportEstimate;
end

for k = 1:0
	sse = SupportEstimate;
	ss = ScreeningLasso(y,X,length(sse));
	SupportEstimate = intersect(SupportEstimate,ss);
	length(SupportEstimate)
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
