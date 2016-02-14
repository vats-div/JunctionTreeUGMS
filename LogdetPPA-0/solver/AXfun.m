%%*********************************************************
%% AXfun: compute AX(k) = <Ak,X>, k = 1:m
%%
%% AX = AXfun(blk,At,X,AL);
%%
%%**********************************************************

  function AX = AXfun(blk,At,X,AL);


  if (nargin < 4); AL = []; end 
  if isempty(AL); existAL = 0; else; existAL = 1; end
%%
  if iscell(At)
     m = size(At{1},2);  
     AX = zeros(m,1);    
     for p = 1:size(blk,1); 
        pblk = blk(p,:);
        if strcmp(pblk{1},'s')
           if (length(pblk{2}) == 1)
              AXtmp = (mexsvec(pblk,X{p})'*At{p,1})';
	   else
              AXtmp = (mexsvec(pblk,sparse(X{p}))'*At{p,1})';
           end
        elseif strcmp(pblk{1},'l') | strcmp(pblk{1},'u')
           AXtmp = (X{p}'*At{p,1})'; 
        end
	AX = AX + AXtmp; 
     end
  else
     if strcmp(blk{1,1},'s')
	if (length(blk{1,2})==1)
           AX = (mexsvec(blk,X)'*At)';
	else
           AX = (mexsvec(blk,sparse(X))'*At)';
        end
     elseif strcmp(blk{1,1},'l') | strcmp(blk{1,1},'u')
        AX = (X'*At)'; 
     end
  end
  if (existAL)
     if strcmp(AL.matfct_options,'chol')
        AX = AL.Rt \ AX; 
     elseif strcmp(AL.matfct_options,'spcholmatlab')
        AX = mexfwsolve(AL.R,AX); 
     end
  end
%%*********************************************************
