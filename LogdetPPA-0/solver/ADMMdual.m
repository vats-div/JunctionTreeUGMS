%%********************************************************************
%% ADMMdual: Alternating direction method of multiplier on the 
%%           dual problem.
%%
%%********************************************************************

  function [X,y,Z,par,iter,sig0] = ADMMdual(blk,At,C,b,L,mu,par,X,y,Z); 

    if ~exist('X');   
       X = ops(blk,'identity'); 
       Z = ops(blk,'identity');
       y = zeros(length(b),1);   
    end
    sig0 = 1; 
    dim = par.dim; 
    if (max(dim(1)) >= 800); 
       maxiter = 20; 
    elseif max(dim(1) >= 400); 
       maxiter = 30; 
    else
       maxiter = 50; 
    end
%%
    AX = AXfun(blk,At,X);   
    AC = AXfun(blk,At,C); 
    AZ = AXfun(blk,At,Z); 
    Rp = b-AX; 
    normX = ops(X,'norm');
    normb = max(1,norm(b)); normC = max(1,ops(C,'norm'));
%%
    for iter = 1:maxiter
       Xold  = X; 
       AXold = AX; 
       AZold = AZ; 
       normXold = normX; 
       par.gamma = sig0.*mu; 
       rhs = Rp/sig0 + AC-AZ; 
       y   = linsysolvefun(L,rhs);
       XmsigC = ops(X,'-',ops(sig0,'*',C)); Aty = Atyfun(blk,At,y); 
       W  = ops(XmsigC,'+',ops(sig0,'*',Aty)); 
       [X,par]  = project3(blk,W,par);
       AX = AXfun(blk,At,X); 
       AZ = (AX-AXfun(blk,At,W))/sig0; 
       Rp = b-AX; 
       normX  = ops(X,'norm'); 
       normZ  = ops(ops(par.phiminus,'/',sig0),'norm');   
       normRp = norm(Rp);
       normRd = sqrt(normX^2+normXold^2-2*blktrace(blk,X,Xold))/sig0; 
       prim_infeas = normRp/normb; 
       dual_infeas = normRd/normC; 
       ratio = prim_infeas/dual_infeas; 
       if (rem(iter,5)==1)
          %fprintf('\n %2.0f  %3.2e  %3.2e %3.2e %3.2e',...
          %        iter,prim_infeas,dual_infeas,ratio,sig0);
       end
       if (rem(iter,5) == 1) %% 2009-Nov-11
          if (iter==1) & (ratio < 0.1 | ratio > 10)
             const = 10; 
             X = Xold; Rp = b-AXold; AZ = AZold;
          else; 
             const = 2; 
          end
          if (ratio < 0.1) 
             sig0 = min(1e3,const*sig0); 
          elseif (ratio > 10) 
             sig0 = max(1e-2,sig0/const); 
          end
       end
       if (max(prim_infeas,dual_infeas) < 1e-2)
          break; 
       end
    end
    Z = ops(ops(X,'-',W),'/',sig0); 
%%********************************************************************
