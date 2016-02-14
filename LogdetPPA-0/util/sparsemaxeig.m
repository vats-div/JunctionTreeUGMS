%%*****************************************************
%% sparsemaxeig: Sparse maximum eigenvalue
%% max {<G,X>  : <I,X> = 1, <E,abs(X)> <= k, X psd}
%%
%%  X = U-V, with U,V >= 0.
%%
%%*****************************************************

   function [blk,At,C,b] = sparsemaxeig(G,k)

   n = length(G);    

   n2 = n*(n+1)/2; 
   blk{1,1} = 's'; blk{1,2} = n; 
   blk{2,1} = 'l'; blk{2,2} = 2*n2+1; 

   II = speye(n2,n2);   
   const = sqrt(2*n2+1);
   ee = svec(blk(1,:),ones(n,n),0);

   At{1,1} = [svec(blk(1,:),speye(n,n),1), II, sparse(n2,1)]; 
   At{2,1} = [sparse(1,2*n2+1); -II, II, sparse(n2,1); ee',ee',1]'; 
   b = [1; zeros(n2,1); k]; 

   C{1,1} = G; 
   C{2,1} = zeros(2*n2+1,1); 
%%*****************************************************
