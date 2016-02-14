%%*****************************************************
%% hardspheres: hard spheres problem
%% min {z : diag(Y) = e, Y_{ij} <= (z+1)-1, Y psd}
%%
%% Optimal Y = (1+a)*eye(n)-a*ones(n), a = 1/(n-1).
%% 
%% Note: z+1 >=0. 
%%*****************************************************

  function [blk,At,C,b] = hardspheres(n)

  blk{1,1} = 's'; blk{1,2} = n;  
  C{1,1} = sparse(n,n); 
  n2 = n*(n+1)/2; 
  Acell = cell(1,n2); 
  for k = 1:n
      Acell{k} = spconvert([k,k,1;n,n,0]); 
  end
  count = n+1; 
  for j = 1:n 
     for i = 1:j-1
        Acell{count} = spconvert([i,j,0.5;j,i,0.5;n,n,0]); 
        count = count + 1; 
     end
  end
  At(1) = svec(blk(1,:),Acell,1);  
  n3 = n*(n-1)/2; 
  b = [ones(n,1); -ones(n3,1)];
   
  blk{2,1} = 'l'; blk{2,2} = n3 + 1; 
  C{2,1} = [1; zeros(n3,1)]; 
  At{2,1} = [sparse(n,n3+1); -ones(n3,1), speye(n3,n3)]'; 
%%*****************************************************

