%%*********************************************************************
%% addcut: Add the constraints
%%
%% G(i,j) >= 0 if G(i,j) < 0.
%%
%%*********************************************************************

  function  [blk2,At2,C2,b2] = addcut(G,blk,At,C,b);

  r2 = sqrt(2);
  n  = length(G); 
%%
  nzidx = [];
  count = 0;
  for j = 2:n
     Gj = G(1:j-1,j);
     idx = find(Gj < -1e-6);
     len = length(idx);
     if (len)
        tmp = [1:len];
        II(count+tmp,1) = idx;
        JJ(count+tmp,1) = j*ones(len,1);
        j2 = (j-1)*j/2;
        tmp1 = j2+idx;
        tmp2 = count+[1:len]';
        nzidx = [nzidx;  [tmp1,tmp2,ones(len,1)]];
        count = count + len;
     end
  end
  Atmp = spconvert([nzidx; n*(n+1)/2, count, 0]);
%%
  blk2 = blk; At2 = At; C2 = C; b2 = b;    
  blk2{2,1} = 'l'; blk2{2,2} = count;
  m = length(b); 
  At2{1,1} = [At{1},Atmp];
  At2{2,1} = [sparse(m,count); -speye(count)]';
  C2{2,1}  = sparse(count,1);
  b2 = [b; zeros(count,1)];
%%*********************************************************************
