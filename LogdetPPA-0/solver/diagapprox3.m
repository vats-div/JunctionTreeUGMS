%%*************************************************************************
%% diagapprox: compute an approximation of 
%% diag (A diag(DPhi) At) 
%%*************************************************************************

   function dd = diagapprox(blk,At,par); 

      m = size(At{1},2); 
      dd = zeros(m,1); 
      for p = 1:size(blk,1);
         pblk = blk(p,:); 
         if strcmp(pblk{1},'s')
            PP = (par.P{p}) .* (par.P{p}); 
            tmp = PP *par.Dsch{p}* PP';
            vv = vectriu(pblk,tmp);
         elseif strcmp(pblk{1},'l')
            vv = par.Dsch{p};
         end
         len = length(vv); 
         D  = spdiags(vv,0,len,len); 
         dd = dd + sum(At{p}.*(D*At{p}))'; 
      end
%%*************************************************************************
%%*************************************************************************
%% vectriu: stack the upper triangular part of a matrix
%%          into a vector. 
%%
%% xvec = vectriu(blk,x);
%%
%%*************************************************************************

  function xvec = vectriu(blk,x); 

  if ~iscell(x)
     r = blk{1,2}; 
     ismt = isempty(x); 
     if strcmp(blk{1,1},'s') & (r > 0);      
        r2 = r*(r+1)/2; 
        xvec = zeros(r2,1); 
        idx = 0; 
        for j = 1:r  
           xvec(idx+[1:j]) = x(1:j,j); 
           idx = idx + j; 
        end
     elseif strcmp(blk{p,1},'l') & (r > 0);
        xvec = x; 
     end
  else
     xvec = [];
     for p = 1:size(blk,1);         
        xvec = [xvec; vectriu(blk(p,:),x{p})]; 
     end   
  end
%%*************************************************************************
