%%**********************************************************************
%% blklogdet: compute logdet(X)
%% 
%%**********************************************************************

  function logdet = blklogdet(blk,X); 

  for p = 1:size(blk,1)
     pblk = blk(p,:);
     n = sum(pblk{2}); 
     if strcmp(pblk{1},'s')
        pert = sqrt(norm(X{p},'fro'))*1e-16; 
        [Xchol,indef] = chol(X{p}+pert*speye(n,n)); 
        logdet(p,1) = 2*sum(log(diag(Xchol)+eps));  
     elseif strcmp(pblk{1},'l')
        logdet(p,1) = sum(log(X{p}+1e-16)); 
     end
  end
%%**********************************************************************




