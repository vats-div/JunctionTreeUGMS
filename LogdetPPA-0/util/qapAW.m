%%******************************************************
%% qapAW: generate SDP data for the problem QAP-AW1 in 
%%        the paper by Povh and Rendl.
%%
%%******************************************************

   function [blk,Ft,CC,bb] = qapAW(A,B,options)
   
   n = size(A,1);
   if (size(B,1)~= n); 
      error('A and B must have the same size'); 
   end
   if (nargin == 2); options = 2; end
   
   n2 = n*n;
   blk{1,1} = 's'; blk{1,2} = n2;
%%
%% 
%%
   CC{1} = kron(B,A); 
   if norm(A-A','fro') | norm(B-B','fro')
      CC{1} = 0.5*(CC{1} + CC{1}'); 
   end
%%
%% sum_i Y^{ii} = I_n
%% create n*(n+1)/2 constraints
%%   
   pp = [0:n:n*(n-1)]'; vv = (1/sqrt(2))*ones(n,1); 
   Fcell = cell(1,n*(n+1)/2);
   count = 0; 
   for j = 1:n
      for i = 1:j
         row = pp + i; col = pp + j; 
         if (i==j)
            Fcell{count+1} = spconvert([row,col,ones(n,1); n2,n2, 0]); 
         else
            Fcell{count+1} = spconvert([row,col,vv; n2,n2, 0]); 
         end
         count = count+1; 
      end
   end   
   Ft = svec(blk,Fcell,1); 
   bb = svec(blk,speye(n,n)); 
test = 0; 
if (test)
   Y = randn(n2,n2); Y = Y+Y';  
   Ysum = zeros(n); 
   for k=1:n; idx = (k-1)*n + [1:n]; Yk = Y(idx,idx); Ysum = Ysum + Yk; end
   blk2{1,1}='s'; blk2{1,2} = n; 
   Ysum2 = smat(blk2,(svec(blk,Y)'*Ft{1})');
   norm(Ysum2./Ysum - ones(n),'fro')
end
%%
%% <I_n,Y^{ij}> = delta_{ij} for j=1:n, i=1:j. 
%% create n*(n+1)/2 constraints
%%   
   pp = [1:n]'; vv = ones(n,1)/2; 
   count = 0; 
   btmp = zeros(n*(n+1)/2,1); 
   for j = 1:n
      for i = 1:j
         row = pp + (i-1)*n; col = pp + (j-1)*n; 
         if (i==j) 
            Fcell{count+1} = spconvert([row,col,ones(n,1); n2,n2, 0]); 
            btmp(count+1) = 1; 
            count = count+1; 
         elseif (i < j)
            Fcell{count+1} = spconvert([row,col,vv; n2,n2, 0]); 
            btmp(count+1) = 0; 
            count = count+1; 
         end
      end
   end
   if (options == 1)
      %%
      %% add <E,Y> = n2
      %%
      idx = [1:n*(n+1)/2-1]; 
      alpha = n2; 
      ee = svec(blk,ones(n2,n2),1); 
      Ftmp = svec(blk,Fcell(idx),1); 
      Ft{1} = [Ft{1}, Ftmp{1}, ee/alpha];
      bb = [bb; btmp(idx); n2/alpha];  
   else
      %%
      %% add <ee^T,Y^{ij}> = 1 for j=1:n, i=1:j. 
      %% 
      Fcell2 = cell(1,n*(n+1)/2);
      alpha = sqrt(n);
      btmp2 = ones(n*(n+1)/2,1)/alpha;      
      rr = [1:n]'*ones(1,n); rr = rr(:);  
      cc = ones(n,1)*[1:n];  cc = cc(:);
      vv = ones(n2,1)/(2*alpha);
      ee = ones(n2,1)/alpha;
      count = 0; 
      for j = 1:n
         for i = 1:j
	    row = rr + (i-1)*n; col = cc + (j-1)*n;
	    if (i==j) 
               Fcell2{count+1} = spconvert([row,col,ee;n2,n2,0]); 
	       count = count+1; 
	    elseif (i < j)
               Fcell2{count+1} = spconvert([row,col,vv;n2,n2,0]); 
   	       count = count+1; 
            end
         end
      end
      idx = [1:n*(n+1)/2-1];
      Ftmp = svec(blk,Fcell(idx),1); Ftmp2 = svec(blk,Fcell2(idx),1); 
      Ft{1} = [Ft{1}, Ftmp{1}, Ftmp2{1}];
      bb = [bb; btmp(idx); btmp2(idx)];  
   end


if (test)
   for k = 1:n;
      for j = 1:k;
         idxj = (j-1)*n + [1:n]; 
         idxk = (k-1)*n + [1:n]; 
         Yjk = Y(idxj,idxk); trY(j,k) = trace(Yjk);    
         trY2(j,k) = svec(blk,Y)'*Ftmp{1}(:,j+k*(k-1)/2);
      end
   end
   norm(triu(trY-trY2),'fro')
end
if (test & options==2)
   for k = 1:n;
      for j = 1:k;
         idxj = (j-1)*n + [1:n]; 
         idxk = (k-1)*n + [1:n]; 
         Yjk = Y(idxj,idxk); Z(j,k) = sum(sum(Yjk));    
         Z2(j,k) = svec(blk,Y)'*Ftmp2{1}(:,j+k*(k-1)/2);
      end
   end
   norm(triu(Z-Z2),'fro')
end

%%******************************************************
