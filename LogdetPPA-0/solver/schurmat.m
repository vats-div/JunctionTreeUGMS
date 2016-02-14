%%***************************************************************************
%% schurmat: 
%%
%%***************************************************************************

  function MM = schurmat(blk,At,W,gamma); 

  if (nargin <= 3); error('schurmat: need at least 3 inputs'); end
  if (nargin < 4); gamma = 0; end
  
  m = size(At{1},2);
%%
  par.gamma = gamma; 
  [X,par]  = project3(blk,W,par);
%%
  MM = zeros(m,m);  
  for p = 1:size(blk,1)
     pblk = blk(p,:); 
     if strcmp(pblk{1},'s')
        P = par.P{p};  Pt = P'; 
        Mp = zeros(m,m); 
        for k = 1:m
           Ak = smat(pblk,At{1}(:,k)); 
           tmp = Pt*Ak*P; tmp = 0.5*(tmp+tmp'); 
           B = P*(par.Dsch{p}.*tmp)*Pt; 
           Mp(:,k) = (svec(pblk,B)'*At{1})';
        end
     elseif strcmp(pblk{1},'l')
        n = sum(pblk{2}); 
        Mp = At{p}'*spdiags(par.Dsch{p},0,n,n)*At{p}; 
     end
     MM = MM + Mp; 
  end 
  MM = 0.5*(MM+MM'); 
%%***************************************************************************
