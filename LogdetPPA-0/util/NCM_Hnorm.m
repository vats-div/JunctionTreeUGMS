%%*****************************************************
%% NCM: nearest correlation matrix
%% min { sum(W.*abs(X-G)) : diag(X) = e, X psd}
%%
%% X-G = U-V, with U,V >= 0.
%%
%% W = a nonnegative weight matrix. 
%%*****************************************************

   function [blk,At,C,b] = NCM_Hnorm(G,W)

   n = length(G); 
   if (nargin == 1); W = ones(n); end

   n2 = n*(n+1)/2; 
   blk{1,1} = 's'; blk{1,2} = n; 
   blk{2,1} = 'l'; blk{2,2} = 2*n2; 

   ee = svec(blk(1,:),ones(n,n)); 
   ww = svec(blk(1,:),W); 
   II = speye(n2,n2); 
   gg = svec(blk(1,:),G); 

   Acell = cell(1,n); 
   for k = 1:n; Acell{k} = spconvert([k,k,1;n,n,0]); end
   Atmp = svec(blk(1,:),Acell,1); 

   At{1,1} = [Atmp{1}, II]; 
   At{2,1} = [sparse(n,2*n2); -II, II]'; 
   b = [ones(n,1); gg]; 

   C{1,1} = sparse(n,n); 
   C{2,1} = [ww; ww]; 
%%*****************************************************
