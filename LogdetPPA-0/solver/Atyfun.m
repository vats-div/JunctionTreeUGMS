%%*********************************************************
%% Atyfun: compute sum_{k=1}^m yk*Ak. 
%%
%%  Q = Atyfun(blk,At,y,AL);
%%**********************************************************

  function Q = Atyfun(blk,At,y,AL);

  if (nargin < 4); AL = []; end
  if isempty(AL); existAL = 0; else; existAL = 1; end; 

  isspAy = ones(size(blk,1),1); 
  Q = cell(size(blk,1),1);
%%
  if (existAL)
     if strcmp(AL.matfct_options,'chol')
        y = AL.R \ y;
     elseif strcmp(AL.matfct_options,'spcholmatlab')
        y = mexbwsolve(AL.Rt,y); 
     end
  end
  if iscell(At)
     for p = 1:size(blk,1)
        pblk = blk(p,:);
        if strcmp(pblk{1},'s')
           Q{p} = mexsmat(pblk,At{p,1}*y);
        elseif strcmp(pblk{1},'l') | strcmp(pblk{1},'u')
           Q{p} = At{p,1}*y; 
        end 
     end
  else
     if strcmp(blk{1,1},'s')
        Q = mexsmat(blk,At*y);
     elseif strcmp(blk{1,1},'l') | strcmp(blk{1,1},'u')
        Q = At*y;
     end
  end
%%********************************************************* 

