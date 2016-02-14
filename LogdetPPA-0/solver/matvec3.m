%%*************************************************************************
%% matvec: matrix-vector multiply.
%% matrix B = sig* (A*DPi*At) = sig * A(PxP) Dsch (PtxPt)At 
%%*************************************************************************

   function By = matvec(blk,At,par,y);
   
   N = length(y); 
   if (norm(y) == 0); By = zeros(N,1); return; end
 %%
   By = zeros(N,1);
   for p = 1:size(blk,1)
      pblk = blk(p,:);
      Aty  = Atyfun(pblk,At{p},y); 
      if strcmp(pblk{1},'s')
         tmp1 = par.Pt{p}*Aty*par.P{p};
         tmp2 = par.Dsch{p}.*tmp1;
         tmp1 = par.P{p}*tmp2*par.Pt{p}; 
         tmp2 = 0.5*(tmp1 + tmp1');
         By = By + par.sig*AXfun(pblk,At{p},tmp2); 
      elseif strcmp(pblk{1},'l')
         tmp = par.Dsch{p}.*Aty; 
 	 By = By + par.sig*(tmp'*At{p})';
      end
   end
%%
   if (par.use_proximal); 
      By = By + (par.H2.*y)/par.sighat; 
   else
      sighat = max([1e4,10*par.sig]); 
      By = By + y/sighat; 
   end
   if isfield(par,'epsilon')
      By = By + par.epsilon*y; 
   end
%%*************************************************************************



