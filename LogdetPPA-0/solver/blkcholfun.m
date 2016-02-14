%%******************************************************************
%% blkcholfun: compute Cholesky factorization of X. 
%%          
%%  [Xchol,indef] = blkcholfun(blk,X,permX); 
%%  
%%  X = Xchol'*Xchol;
%%
%%******************************************************************
 
  function [Xchol,indef] = blkcholfun(blk,X); 

  if iscell(X) 
     Xchol = cell(size(X));
     indef = zeros(size(blk,1),1);   
     for p = 1:size(blk,1) 
        pblk = blk(p,:); 
        if strcmp(pblk{1},'s');
           [Xchol{p},indef(p)] = chol(X{p}); 
        elseif strcmp(pblk{1},'l'); 
           if any(X{p} <= 0) 
              indef(p) = 1;
           else
              Xchol{p} = sqrt(X{p}); 
           end 
        end
     end
  else
     if strcmp(blk{1},'s')
        [Xchol,indef] = chol(X); 
     elseif strcmp(blk{1},'l')              
        if (any(X) <= 0)
           Xchol = []; 
           indef = 1; 
        else
           Xchol = sqrt(X); 
           indef = 0; 
        end
     end
  end 
%%******************************************************************     


