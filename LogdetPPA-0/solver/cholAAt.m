%%***********************************************************
%% cholAAt: compute Cholesky factorization of A*At. 
%%
%%***********************************************************

   function [L,AAt] = cholAAt(blk,At,m); 

   AAt = sparse(m,m);  
   for p = 1:size(blk,1)
      AAt = AAt + At{p,1}'*At{p,1}; 
   end
   pertdiag = 1e-13*ones(m,1); 
   AAt = AAt + spdiags(pertdiag,0,m,m);
%%
   if (nnz(AAt) < 0.2*m*m); use_spchol=1; else; use_spchol=0; end
%%  
   if (use_spchol)
      [L.R,L.p,L.perm] = chol(AAt,'vector'); 
      L.Rt = L.R'; 
      L.matfct_options = 'spcholmatlab'; 
   else
      if issparse(AAt); AAt = full(AAt); end;           
      L.matfct_options = 'chol';
      L.perm = [1:m]; 
      [L.R,indef] = chol(AAt); 
   end
%%***********************************************************
