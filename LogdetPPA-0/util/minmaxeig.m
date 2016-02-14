%%*********************************************************************
%% minmaxeig: 
%%
%% minimize maxeig(B_0 + x_1*B_1 + ... + x_m*B_m)
%%
%%
%% The problem is equivalent to the following SDP:
%% 
%%    minimize{ t : t*I - B(x) >= 0 }
%%
%% with B(x) = B_0 + x_1*B_1 + ... + x_m*B_m.
%%
%% Adapted from Vandenberghe and Boyd, SIAM Review 96.        
%%-------------------------------------------------------------------- 
%%
%%  [blk,At,C,b] = minmaxeig(n,m); 
%%
%%*********************************************************************

   function [blk,At,C,b] = minmaxeig(n,m);   

   tmp = randn(n); B{1} = tmp*tmp';
   for k = 1:m; tmp=randn(n); B{k+1}=tmp+tmp'; end
%%
%% A_i = B_i, i=1,...,m; A_{m+1} = I
%%     
   blk{1,1} = 's'; blk{1,2} = n;
   b = [zeros(m,1); -1];
   C = -B{1}; 
   A = cell(1,m+1); 
   for k = 1:m; A{k} = B{k+1}; end
   A{m+1} = -speye(n,n); 
   At = svec(blk,A,1); 
%%*********************************************************************
