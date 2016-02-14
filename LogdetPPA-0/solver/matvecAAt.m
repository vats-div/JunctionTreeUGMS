%%************************************************************************
%% matvecAAt: matrix-vector multiply.
%% matrix H = A*At
%%************************************************************************

  function Hy = matvecAAt(blk,At,par,y,AL)

  if (nargin < 5); AL = []; end

  N = length(y);
  if (norm(y) == 0); Hy = zeros(N,1); return; end

  tmp = Atyfun(blk,At,y,AL); 
  Hy = AXfun(blk,At,tmp,AL);
%%************************************************************************
