%%****************************************************************************************
%% logdetPPA: The software is designed to solve standard primal and 
%%            dual log-determinant SDPs of the form:
%%
%% (P)  min  (C1,X1) + ... + (CN,XN) - mu1*logdet(X1) - ... - muN*logdet(XN) 
%%      s.t. A1(X1) + ... + AN(XN) = b, 
%%           X1,...,XN positive definite
%%
%% (D)  max  (b,y) + mu1*logdet(Z1) + ... + muN*logdet(ZN) + constant 
%%      s.t. A1t(y) + Z1 = C1, ... ,ANt(y) + ZN = CN, 
%%           Z1,..,ZN positive definite
%%
%%      where Xk, Zk are either symmetric positive definite matrices or positive vectors.  
%%      The parameters mu1,...,muN are given positive scalars. 
%%      The constant term in (D) depends on the demensions of X1,...,XN and mu1,...,muN. 
%%
%% [obj,X,y,Z,info,runhist] = logdetPPA(blk,At,C,b,mu,OPTIONS);
%%
%% Input: blk,At,C,b (problem data whose input format is the same as in SDPT3) 
%%        mu = [mu1,mu2,...,muN]; 
%% Output: obj  = [primal objective value, dual objective value]
%%         X    = a cell array containing the primal variable.
%%         y    = a vector of Lagrange multiplers.
%%         Z    = a cell array containing the dual slack variable.
%%         info = a structure array containing various information such as,
%%                iter    : number of outer iteration
%%                itersub : total number of sub-problems solved 
%%                obj     : objective value 
%%                relgap  : relative duality gap
%%                pinfeas : primal infeasibility 
%%                dinfeas : dual infeasibility 
%%                termcode: termination code 
%%                          (0 : problem is solved to required tolerance, tol)
%%                          (-1: problem is partially solved to the required tolerance
%%                               with pinfeas < sqrt(tol) and dinfeas < tol)
%%                          ( 1: problem is not solved successful due to stagnation) 
%%                          ( 2: maximum number of iterations reached)
%%        runhist = a structure array containing the history of the run. 
%%****************************************************************************************
%% LogdetPPA: version 0
%% Copyright (c) 2009 by
%% Chengjing Wang, Defeng Sun, and Kim-Chuan Toh
%% Last modified: 27 May 1010
%%****************************************************************************************
 
   function [obj,X,y,Z,info,runhist] = logdetPPA_newton(blk,At,C,b,mu0,OPTIONS);

   if (nargin < 6); OPTIONS = []; end

   tol = 1e-6; 
   sig = 10;   
   maxiter    = 100; 
   maxitersub = 30;
   precond    = 1; 
   maxitpsqmr = 100; 
   stagnate_check_psqmr = 0; 
   scale_data = 2; 
   printlevel = 1; 
   plotyes    = 0; 
   if (norm(mu0) > 1e-12)
      use_smoothing = 1; 
   else
      use_smoothing = 0; 
   end
   use_proximal = 1; %% better choice: 1;      
   switch_alt_newton_tol = 1e-2; 
%%
   if isfield(OPTIONS,'tol');        tol           = OPTIONS.tol; end
   if isfield(OPTIONS,'sigma');      sig           = OPTIONS.sigma; end
   if isfield(OPTIONS,'maxiter');    maxiter       = OPTIONS.maxiter; end
   if isfield(OPTIONS,'maxitersub'); maxitersub    = OPTIONS.maxitersub; end   
   if isfield(OPTIONS,'precond');    precond       = OPTIONS.precond; end   
   if isfield(OPTIONS,'maxitpsqmr'); maxitpsqmr    = OPTIONS.maxitpsqmr; end   
   if isfield(OPTIONS,'stagnate_check_psqmr'); 
      stagnate_check_psqmr = OPTIONS.stagnate_check_psqmr; 
   end   
   if isfield(OPTIONS,'scale_data'); scale_data    = OPTIONS.scale_data; end   
   if isfield(OPTIONS,'printlevel'); printlevel    = OPTIONS.printlevel; end   
   if isfield(OPTIONS,'plotyes');    plotyes       = OPTIONS.plotyes; end  
   if isfield(OPTIONS,'proximal');   use_proximal  = OPTIONS.proximal; end 
   if isfield(OPTIONS,'smoothing');  use_smoothing = OPTIONS.smoothing; end
   if isfield(OPTIONS,'switch_alt_newton_tol'); 
      switch_alt_newton_tol = OPTIONS.switch_alt_newton_tol; 
   end
%%
   if ~iscell(C); tmp = C; clear C; C{1} = tmp; end
%%
   dim = zeros(1,2); numblk = zeros(1,2); 
   sdpblkidx = [];    
   for p = 1:size(blk,1)
      pblk = blk(p,:); 
      if strcmp(pblk{1},'s')
         dim(1) = dim(1) + sum(pblk{2});
         numblk(1) = numblk(1) + length(pblk{2}); 
         sdpblkidx = [sdpblkidx, p]; 
      elseif strcmp(pblk{1},'l')
         dim(2) = dim(2) + sum(pblk{2});
         numblk(2) = numblk(2) + length(pblk{2}); 
      end
   end   
   sdpblkidx = [sdpblkidx, setdiff([1:size(blk,1)],sdpblkidx)];
   blk = blk(sdpblkidx,:); At = At(sdpblkidx); C = C(sdpblkidx); 
   if (printlevel)
      fprintf('\n num. of constraints = %2.0d',length(b));      
      if dim(1); 
         fprintf('\n dim. of sdp  var = %2.0d,',dim(1)); 
         fprintf('   num. of sdp  blk = %2.0d',numblk(1)); 
      end
      if dim(2); fprintf('\n dim. of linear var = %2.0d',dim(2)); end
   end
   for p = 1:size(blk,1); nblk(p,1) = sum(blk{p,2}); end
   if size(mu0,1) < size(mu0,2); mu0 = mu0'; end
   n = ops(C,'getM');  
   m = length(b); 
   tstart = clock; 
%%
%% scaling
%%
   b0 = b; C0 = C; 
   if (scale_data)
      normA0 = min(1e12,max(1,ops(At,'norm'))); 
      normb0 = max(1,norm(b)); normC0 = min(1e12,max(1,ops(C,'norm'))); 
      for p = 1:size(blk,1)
         At{p} = At{p}/normA0; 
      end
      C = ops(C,'/',normC0); 
      b = b/normb0; 
      objscale = normb0*normC0/normA0; 
   else
      normA0 = 1; normb0 = 1; normC0 = 1;
      objscale = 1;  
   end
   normb = max(1,norm(b)); normC = max(1,ops(C,'norm'));
   if (scale_data)
      mu = mu0*normA0/(normb0*normC0);
   else
      mu = mu0; 
   end
%%
%% use an alternating direction method to compute 
%% initial X,y,Z  
%%
   X = ops(blk,'identity'); 
   Z = ops(blk,'identity');
   y = zeros(m,1);   
   [L,AAt] = cholAAt(blk,At,m); 
   normX = ops(X,'norm');
   normZ = ops(Z,'norm');  
   AX = AXfun(blk,At,X);   
   Rp = b-AX; 
   if (max(dim(1)) >= 300); 
      maxiterAlt = 20; 
   elseif max(dim(1) >= 200); 
      maxiterAlt = 30; 
   else
      maxiterAlt = 50; 
   end
   sig0 = 1;
   iter = 0;
   if (scale_data==2)
      AC = AXfun(blk,At,C); 
      AZ = AXfun(blk,At,Z); 
      for iter = 1:maxiterAlt
         Xold  = X; 
         AXold = AX; 
         AZold = AZ; 
         normXold = normX; 
         par.gamma = sig0.*mu; 
         rhs = Rp/sig0 + AC-AZ; 
         y   = linsysolvefun(L,rhs);
         XmsigC = ops(X,'-',ops(sig0,'*',C)); Aty = Atyfun(blk,At,y); 
         W     = ops(XmsigC,'+',ops(sig0,'*',Aty)); 
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
         if rem(iter,5)==1
            fprintf('\n %2.0f  %3.2e  %3.2e %3.2e %3.2e',...
                    iter,prim_infeas,dual_infeas,ratio,sig0);
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
   end
%%
   if (scale_data==2) 
      if (normZ > 1e-2)
         norm_ratio = full(normX/normZ);
      else
         normtmp = max(1e-2,ops(ops(C,'-',Aty),'norm')); 
         norm_ratio = full(normX/normtmp);
      end
      X = ops(X,'/',norm_ratio); y = y/norm_ratio;
      for p = 1:size(blk,1)
         At{p} = At{p}*norm_ratio; 
      end
      normA0 = normA0/norm_ratio;
      objscale = objscale*norm_ratio;       
      mexScalarxMatrix(AAt,norm_ratio^2);
      mexScalarxMatrix(L.R,norm_ratio); 
      if isfield(L,'Rt'); 
         mexScalarxMatrix(L.Rt,norm_ratio); 
      end
   end
%% 
   if (scale_data==2) 
      sig = 10;
   else
      sig = max(sig,max(1e-2,normX)/max(1e-2,normZ)); 
   end
   if (scale_data)
      mu = mu0*normA0/(normb0*normC0);
   else
      mu = mu0; 
   end
   mutarget = mu; 
%%
   XmsigC = ops(X,'-',ops(sig,'*',C)); Aty = Atyfun(blk,At,y);
   W   = ops(XmsigC,'+',ops(sig,'*',Aty)); 
   Xold = X; 
   normXold = normX; 
   par.gamma = sig.*mu; 
   [X,par]  = project3(blk,W,par);
   normX    = ops(X,'norm'); 
   AX   = AXfun(blk,At,X);
   Rp   = b-AX; 
   normRp = norm(b-AX); 
   normRd = sqrt(normX^2+normXold^2-2*blktrace(blk,X,Xold))/sig; 
   if (use_smoothing)
      pobjadd = -sum(nblk.*mu*log(normb0/normA0));
      dobjadd = sum(nblk.*mu.*log(normC0)) + sum(nblk.*mu.*(1-log(mu0+eps)));
      logdetX = blklogdet(blk,X);
      logdetZ = ops(ops(par.phiminus,'/',sig),'sumlog');
      obj(1)  = objscale*(blktrace(blk,C,X) - mu'*logdetX + pobjadd);
      obj(2)  = objscale*(b'*y + mu'*logdetZ + dobjadd);
   else
      obj = objscale*[blktrace(blk,C,X), b'*y];  
   end
   gap = obj(1) - obj(2);
   prim_infeas = normRp/normb;       
   dual_infeas = normRd/normC; 
   ttime(1) = etime(clock,tstart);
   runhist.prim_obj(1)    = obj(1);
   runhist.dual_obj(1)    = obj(2);
   runhist.gap(1)         = gap; 
   runhist.prim_infeas(1) = prim_infeas; 
   runhist.dual_infeas(1) = dual_infeas; 
   runhist.cputime(1) = ttime; 
   runhist.sigma(1)   = sig; 
   runhist.itersub(1) = iter; 
%%
   fprintf('\n*****************************************************')
   fprintf('********************************')
   fprintf('\n  logdetPPA: a Newton-CG Primal Proximal-Point Algorithm')
   fprintf('\n  use_alternating_newton = 1, use_smoothing = %1.1d, use_proximal = %1.1d',...
                use_smoothing,use_proximal); 
   fprintf('\n  scale_data = %1.1d, precond = %1.1d',scale_data,precond); 
   fprintf('\n  normC = %2.1e, normA = %2.1e, normb = %2.1e',...
           ops(C,'norm'),ops(At,'norm'),norm(b));
   fprintf('\n  number of initial iterations = %2.0f',iter); 
   fprintf('\n  sig0 = %3.2e,  iter = %2.1d',sig0,iter); 
   if (scale_data)
      fprintf('\n  orig-normX = %2.1e, orig-normZ = %2.1e',normX,normZ);  
      fprintf('\n  normX = %2.1e, normZ = %2.1e',ops(X,'norm'),ops(Z,'norm'));
   end
   fprintf('\n ****************************************************'); 
   fprintf('********************************')
   fprintf('\n  it   pinf    dinf          pobj           dobj')
   fprintf('      gap      time   sigma  sigmahat')
   fprintf('\n*****************************************************')
   fprintf('********************************')
   fprintf('\n %3.1d  %3.1e %3.1e  %- 8.7e %- 8.7e %- 3.2e',...
           0,prim_infeas,dual_infeas,obj(1),obj(2),gap); 
   fprintf(' %4.1f  %3.1e',ttime,sig);  
   if (max(prim_infeas,dual_infeas) < tol)
      fprintf('\n--------------------------------------------------------') 
      fprintf('\n max(prim_infeas,dual_infeas) < %3.2e',tol); 
      fprintf('\n--------------------------------------------------------')
      return; 
   end
%%
%% main 
%%
   msg = '';
   diagAAt = max(1e-4,full(diag(AAt))); 
   sighat = 0; 
   if (use_proximal); 
      sighatmin = 1e5; sighat = sighatmin; sighatfac = 10; 
      H = ones(m,1); 
   end
   yhat = y; 
   RpGradratio = 1; 
   par.b = b;
   clear AAt L; 
   if isfield(par,'P'); par = rmfield(par,'P'); end
   if isfield(par,'Dsch'); par = rmfield(par,'Dsch'); end
%%
   for iter = 1:maxiter 
      par.sig = sig;
      par.use_proximal  = use_proximal;
      par.use_smoothing = use_smoothing; 
      if (use_proximal)
         yhat = y;
         sighatmin = max(1,RpGradratio)*sighatmin; %% add: 27-Mar-2008 
         sighatfac = max(1,RpGradratio^2)*sighatfac; 
         sighat = max(sighatmin,sighatfac*sig);    %% add: 07-Apr-2008 
         H = 1./max(1,sqrt(abs(yhat)));  %% better than H = ones(m,1); 
         H2 = H.*H; 
         par.yhat = yhat; par.sighat = sighat; par.H = H; par.H2 = H2; 
      end
      if (use_smoothing)
         mu = max(mutarget,0.4*mu);  %% added: 31-Jul-2009
         par.mu = mu; 
         par.gamma   = par.sig.*par.mu; 
         par.project = 'project3'; 
         matvecfname = 'matvec3'; 
         diagapproxfname = 'diagapprox3';
      else
         par.project = 'project'; 
         matvecfname = 'matvec'; 
         diagapproxfname = 'diagapprox';
      end
      if (dual_infeas < 1e-5) 
         maxitersub = max(maxitersub,50); 
      elseif (dual_infeas < 1e-3)
         maxitersub = max(maxitersub,40); 
      elseif (dual_infeas < 1e-1) 
         maxitersub = max(maxitersub,30); 
      end 
      %%-----------------------------------
      %% Newton method for inner subproblem
      %%-----------------------------------
      if (prim_infeas < 0.2*dual_infeas)
         const = 0.5; 
      elseif (prim_infeas > dual_infeas) 
         const = 0.1;
      else
         const = 0.2;  
      end      
      tolsub   = min(1,const*normRd);  
      printsub = 1; 
      breakyes = 0; 
      %%-----------------------------------
      %% new X, compute W
      %%-----------------------------------
      Xold = X;  
      normXold = normX;
      XmsigC = ops(X,'-',ops(sig,'*',C)); Aty = Atyfun(blk,At,y); 
      W   = ops(XmsigC,'+',ops(sig,'*',Aty));       
      [X,par] = feval(par.project,blk,W,par); 
      ysub = y; 
      clear subhist
      subhist.psqmr(1)    = 0; 
      subhist.findstep(1) = 0; 
      subhist.y(:,1)      = ysub;  
      %%-----------------------------------------
      %% inner iteration
      %%-----------------------------------------
      for itersub = 1:maxitersub

         normX = ops(X,'norm'); 
         Ly = b'*ysub - normX^2/(2*sig); 
         AX = AXfun(blk,At,X);
         Rpsub  = b-AX;
	 GradLy = Rpsub; 
         normRpsub = norm(Rpsub);  
	 if (use_proximal)
	    Ly = Ly - norm(H.*(ysub-yhat))^2/(2*sighat);
            GradLy = GradLy - H2.*(ysub-yhat)/sighat;
         end
	 if (use_smoothing)
            Ly = Ly - sum(par.mu.*par.logphi); 
         end         
         RpGradratio = max(1,norm(Rpsub)/norm(GradLy)); 
         %%-----------------------------------------
         normRdsub = sqrt(normX^2+normXold^2-2*blktrace(blk,X,Xold))/sig;
	 const = 0.2; 
         tolsub = max(min(1,const*normRdsub),1e-2*tol);
         %%
	 priminf_sub = normRpsub/normb; 
         dualinf_sub = normRdsub/normC; 
         subhist.priminf(itersub) = priminf_sub; 
         subhist.dualinf(itersub) = dualinf_sub; 
         subhist.Ly(itersub)      = Ly; 
	 if (printsub)
            fprintf('\n      %2.0d  %- 11.10e %3.2e %3.2e %3.1f',...
            itersub,Ly,priminf_sub,dualinf_sub,const);
         end
         %%-----------------------------------------
	 if (norm(GradLy) < tolsub) & (itersub > 1) 
	    msg = 'good termination:'; 
	    fprintf('\n       %s  ',msg); 
	    fprintf(' normRp=%3.2e, gradLy=%3.2e, tolsub=%3.2e',...
               normRpsub,norm(GradLy),tolsub);
            breakyes = -1; 
            break; 
         end
         if (max(priminf_sub,dualinf_sub) < tol) 
            msg = sprintf('max(prim_infeas,dual_infeas) < %3.2e',tol);  
	    fprintf('\n       %s',msg);
            breakyes = -2; 
	    break;
         end
         %%----------------------------------------- 
         %% compute Newton direction
	 %% precond = 1, precondition by diag(A*diag(DPhi)*At)
	 %% precond = 2, precondition by AAt
         %%----------------------------------------- 
	 par.epsilon = min(1e-3,0.1*norm(GradLy)); %% good to add
	 par.precond = precond;
	 if (dual_infeas > 1e-4); par.precond = 1; end
         diag_precond_options = 1;
	 if (diag_precond_options==1)
            diagAAt = max(1e-4,feval(diagapproxfname,blk,At,par)); 
         end
	 if (use_proximal)
	    L.invdiagM = 1./full(sig*diagAAt + H2/sighat); 
	 else
	    L.invdiagM = 1./full(sig*diagAAt); 
         end     
	 if (dual_infeas > 1e-3) | (itersub <= 5)
	    maxitpsqmr = max(maxitpsqmr,200); 
	 elseif (dual_infeas > 1e-4)	 
	    maxitpsqmr = max(maxitpsqmr,300); 
	 elseif (dual_infeas > 1e-5)	 
	    maxitpsqmr = max(maxitpsqmr,400); 
	 elseif (dual_infeas > 5e-6)
	    maxitpsqmr = max(maxitpsqmr,500); 
         end
	 par.minitpsqmr = 3;
	 if (dual_infeas > 1e-4)
	    stagnate_check_psqmr = max(stagnate_check_psqmr,20);
	 else
	    stagnate_check_psqmr = max(stagnate_check_psqmr,50);
         end
	 if (itersub > 3 & all(subhist.solve_ok(itersub-[3:-1]) <= -1)) ...
	    & (dual_infeas < 5e-5) 
	    stagnate_check_psqmr = max(stagnate_check_psqmr,100);
         end
	 par.stagnate_check_psqmr = stagnate_check_psqmr;
         if (itersub > 1) 
	    prim_ratio = priminf_sub/subhist.priminf(itersub-1); 
	    dual_ratio = dualinf_sub/subhist.dualinf(itersub-1); 
	 else
   	    prim_ratio = 0; dual_ratio = 0;              
         end 
         rhs = GradLy;
         if (iter < 2 & itersub < 5) %% added: 15-Jun-2009
	    tolpsqmr = min(1,0.1*norm(rhs)); 
         else
   	    tolpsqmr = min(1e-1,0.1*norm(rhs)); %% old: min(5e-3,0.1*norm(rhs)); 
         end
         if (itersub <= 1)
            const2 = 1; 
         else
            const2 = 0.5; %% better than const2 = 1; 
         end
         if (itersub > 1) & (prim_ratio > 0.5 | priminf_sub > 0.1*subhist.priminf(1))
	    const2 = 0.5*const2;
         end
	 if (dual_ratio > 1.1); const2 = 0.5*const2; end; 
         tolpsqmr = const2*tolpsqmr;  %% smaller tolpsqmr often give better prim_infeas.
         [dy,resnrm,solve_ok] = ...
             psqmr(matvecfname,blk,At,rhs,par,L,tolpsqmr,maxitpsqmr); 
	 psqmriter = length(resnrm)-1; 
	 if (printsub)
             fprintf('| %3.1e %3.1e %3.0d',tolpsqmr,resnrm(end),psqmriter);
	     if (~use_smoothing) & isfield(par,'P1')
                fprintf(' %2.1f %2.0d',const2,size(par.P1{1},2)); %%
             else  
                fprintf(' %2.1f ',const2); 
             end  
         end
         %%----------------------------------------- 
         %% line search 
         %%----------------------------------------- 
         if (itersub <= 3) & (dual_infeas > 1e-4) | (iter <= 3)
	    step_options = 1; 
         else
	    step_options = 2; 
         end
         steptol = 1e-5;
	 [par,ysub,W,X,alpha,iterstep] = ...
              findstep(blk,At,par,ysub,W,X,dy,steptol,step_options); 
	 subhist.solve_ok(itersub) = solve_ok;
	 subhist.psqmr(itersub)    = psqmriter; 
	 subhist.findstep(itersub) = iterstep; 
	 subhist.y(:,itersub+1)    = ysub;  
         Ly_ratio = 1; 
	 if (itersub > 1)
   	    Ly_ratio = (Ly-subhist.Ly(itersub-1))/(abs(Ly)+eps);
         end
	 if (printsub)
            fprintf(' %3.2e %2.0f',alpha,iterstep);
	    if (Ly_ratio < 0); fprintf('-'); end
         end   
         if (plotyes) 
            plotfun(blk,par,resnrm); 
         end
         %%----------------------------------------- 
         %% check for stagnation
         %%----------------------------------------- 
	 if (itersub > 4)
	    idx = [max(1,itersub-3):itersub]; 
	    tmp = subhist.priminf(idx); 
	    ratio = min(tmp)/max(tmp);
	    if (all(subhist.solve_ok(idx) <= -1)) & (ratio > 0.9) ... 
	       & (min(subhist.psqmr(idx)) == max(subhist.psqmr(idx))) ...
               & (max(tmp) < 5*tol)
	       fprintf('#')
	       breakyes = 1; break; 
            end
	    const3 = 0.7; 
   	    priminf_1half  = min(subhist.priminf(1:ceil(itersub*const3))); 
            priminf_2half  = min(subhist.priminf(ceil(itersub*const3)+1:itersub));
	    priminf_best   = min(subhist.priminf(1:itersub-1)); 
	    priminf_ratio  = subhist.priminf(itersub)/subhist.priminf(itersub-1);
            dualinf_ratio  = subhist.dualinf(itersub)/subhist.dualinf(itersub-1);  
	    stagnate_idx   = find(subhist.solve_ok(1:itersub) <= -1);
	    stagnate_count = length(stagnate_idx); 
	    idx2 = [max(1,itersub-7):itersub]; 
            if (itersub >= 10) & all(subhist.solve_ok(idx2) == -1) ...
	       & (priminf_best < 1e-2) & (dualinf_sub < 1e-3)
	       tmp = subhist.priminf(idx2); 
	       ratio = min(tmp)/max(tmp); 
	       if (ratio > 0.5) 
  	          fprintf('##'); breakyes = 2; break; 
               end  
            end
	    if (itersub >= 15) & (priminf_1half < min(2e-3,priminf_2half)) ...
	       & (dualinf_sub < 0.8*subhist.dualinf(1)) & (dualinf_sub < 1e-3) ...
	       & (stagnate_count >= 3) 
  	       fprintf('###'); breakyes = 3; break; 
            end
	    if (itersub >= 15) & (priminf_ratio < 0.1) ...
               & (priminf_sub < 0.8*priminf_1half) ...
	       & (dualinf_sub < min(1e-3,2*priminf_sub)) ...
	       & ((priminf_sub < 2e-3) | (dualinf_sub < 1e-5 & priminf_sub < 5e-3)) ...
	       & (stagnate_count >= 3) 
	       fprintf(' $$')
	       breakyes = 4; break; 
            end
	    if (itersub >=10) & (dualinf_sub > 5*min(subhist.dualinf)) ...
	       & (priminf_sub > 2*min(subhist.priminf)) %% add: 08-Apr-2008
	       fprintf('$$$')
	       breakyes = 5; break;
            end
         end
      end 
      %%----- end of inner iteration -------------------
      %%
      if (itersub == maxitersub) | (breakyes > 0)
	 [priminfmin,idxmin] = min(subhist.priminf); 
	 idx  = find(subhist.priminf < 1.1*priminfmin); 
   	 ysub = subhist.y(:,max(idx));
         XmsigC = ops(Xold,'-',ops(sig,'*',C)); Aty = Atyfun(blk,At,ysub);  
         W    = ops(XmsigC,'+',ops(sig,'*',Aty)); 
         [X,par] = feval(par.project,blk,W,par);
         normX = ops(X,'norm'); 
      end
      y  = ysub;
      AX = AXfun(blk,At,X);
      normRp = norm(b-AX); 
      normRd = sqrt(normX^2+normXold^2-2*blktrace(blk,X,Xold))/sig; 
      if (use_smoothing)
         pobjadd = -sum(nblk.*mu*log(normb0/normA0));
         dobjadd = sum(nblk.*mu.*log(normC0)) +sum(nblk.*mu.*(1-log(mu0+eps)));
         logdetX = blklogdet(blk,X);
         logdetZ = ops(ops(par.phiminus,'/',sig),'sumlog'); 
         obj(1)  = objscale*(blktrace(blk,C,X) - mu'*logdetX + pobjadd);
         obj(2)  = objscale*(b'*y + mu'*logdetZ + dobjadd);
      else
         obj = objscale*[blktrace(blk,C,X), b'*y];  
      end
      gap = obj(1) - obj(2);
      prim_infeas = normRp/normb;       
      dual_infeas = normRd/normC; 
      runhist.prim_obj(iter+1)    = obj(1);
      runhist.dual_obj(iter+1)    = obj(2);
      runhist.gap(iter+1)         = gap; 
      runhist.prim_infeas(iter+1) = prim_infeas; 
      runhist.dual_infeas(iter+1) = dual_infeas; 
      runhist.cputime(iter+1) = etime(clock,tstart);
      runhist.sigma(iter+1)   = sig; 
      runhist.itersub(iter+1) = itersub-1;  
      runhist.rank(iter+1)    = length(par.posidx{1}); 
      runhist.psqmr(iter+1)   = sum(subhist.psqmr); 
      runhist.findstep(iter+1)= sum(subhist.findstep); 
      %%
      fprintf('\n %3.0d  %3.1e %3.1e  %- 8.7e %- 8.7e %- 3.2e',...
          iter,prim_infeas,dual_infeas,obj(1),obj(2),gap); 
      fprintf(' %4.1f  %3.1e %3.1e ',runhist.cputime(iter+1),sig,sighat); 
      %%----------------------------------
      %% check for termination
      %%----------------------------------
      idx = [max(iter-1,1):iter]; 
      prim_infeas_ratio = runhist.prim_infeas(idx+1)./runhist.prim_infeas(idx);
      dual_infeas_ratio = runhist.dual_infeas(iter+1)/runhist.dual_infeas(iter);
      if (max(prim_infeas,dual_infeas) < tol)
         msg = sprintf('max(prim_infeas,dual_infeas) < %3.2e',tol); 
         break;       
      elseif ((dual_infeas < 0.5*tol) & all(prim_infeas_ratio > 0.7)) ... 
         | ((prim_infeas_ratio(end) > 0.7) & all(runhist.dual_infeas(iter:iter+1) < 0.5*tol))
         msg = sprintf('lack of progress in prim_infeas'); 
         break;
      elseif (dual_infeas < 0.5*tol) & all(dual_infeas_ratio > 0.7) 
         msg = sprintf('dual_infeas deteriorated'); 
         break;
      elseif (dual_infeas < 0.8*tol) & (breakyes)
         msg = sprintf('dual_infeas < %3.2e',0.8*tol); 
         break; 
      end      
      %%----------------------------------
      %% outer Newton iteration
      %%----------------------------------      
      use_alt_newton = (max(prim_infeas,dual_infeas) < switch_alt_newton_tol); 
      if (use_alt_newton) & (breakyes ~= -2)
         fprintf('\n-----use alternating Newton: ')
         RHS = ops(X,'-',Xold);    
         tolnewton = 5e-2*ops(RHS,'norm'); 
         [dX,ddyy,resnrm2,num_innerPCG] = ...
              psqmr_newton(matvecfname,blk,At,RHS,par,L,tolnewton,5);
         X = ops(Xold,'+',dX); 
         normX = ops(X,'norm'); 
         fprintf('numiter = %1.0f, res_ratio = %3.2e, num_innerPCG = %2.0f',...
         length(resnrm2)-1,resnrm2(end)/resnrm2(1),sum(num_innerPCG));
         if (dual_infeas < 1e-6); y = y + ddyy; end
      end 
      %%----------------------------------
      %% modify sig
      %%----------------------------------
      if (dual_infeas > 1e-5)
         sigtol = 0.5; 
      else
         sigtol = 0.5; 
      end
      if (dual_infeas_ratio > sigtol) 
         const4 = 2.0; 
         sig = min(1e8,const4*sig);
      end
   end   
%%
   if (iter == maxiter); breakyes = 6; msg = 'maximum iteration reached'; end
   Z = ops(ops(X,'-',W),'/',sig);   
   normXscaled = ops(X,'norm'); 
   normyscaled = norm(y); 
   normZscaled = ops(Z,'norm'); 
   if (scale_data)
      X = ops(X,'*',normb0/normA0); y = y*(normC0/normA0);      
      Z = ops(normC0,'*',Z);  
      AX = AXfun(blk,At,X); 
      Rp = b0-normA0*AX; 
      Rd = ops(C0,'-',ops(Atyfun(blk,At,normA0*y),'+',Z));    
      prim_infeas = norm(Rp)/max(1,norm(b0));
      dual_infeas = ops(Rd,'norm')/max(1,ops(C0,'norm'));        
   end
   rel_gap = (obj(1)-obj(2))/max(1,mean(abs(obj)));
   info.iter      = iter;
   info.itersub   = sum(runhist.itersub(2:end)); 
   info.ave_psqmr = sum(runhist.psqmr(2:end))/info.itersub; 
   info.obj     = [obj(1),obj(2)]; 
   info.gap     = obj(1)-obj(2); 
   info.relgap  = rel_gap; 
   info.pinfeas = prim_infeas;
   info.dinfeas = dual_infeas; 
   info.normX   = ops(X,'norm');
   info.normy   = norm(y); 
   info.normZ   = ops(Z,'norm'); 
   info.msg     = msg; 
   if (max(prim_infeas,dual_infeas) < tol)   
      info.termcode = 0; 
   elseif (prim_infeas < sqrt(tol) & dual_infeas < tol); 
      info.termcode = -1;
   else
      if (breakyes > 0 & breakyes < 6);  info.termcode = 1; end
      if (breakyes == 6); info.termcode = 2; end 
   end
   fprintf('\n--------------------------------------------------------') 
   fprintf('------------------------')
   fprintf('\n %s',msg); 
   fprintf('\n--------------------------------------------------------')
   fprintf('------------------------')
   fprintf('\n primal objval = %9.8e',obj(1));
   fprintf('\n dual   objval = %9.8e',obj(2));
   fprintf('\n gap           = %3.2e',gap);
   fprintf('\n relative gap  = %3.2e',rel_gap);
   fprintf('\n prim_infeas   = %3.2e',prim_infeas);
   fprintf('\n dual_infeas   = %3.2e',dual_infeas);
   fprintf('\n CPU time      = %3.2f',runhist.cputime(end)); 
   fprintf('\n number of subproblems = %3.0f',sum(runhist.itersub(2:end))); 
   fprintf('\n average number of psqmr step per subproblem = %3.2f',...
           sum(runhist.psqmr(2:end))/sum(runhist.itersub(2:end))); 
   fprintf('\n average number of linesearch per subproblem = %3.2f',...
           sum(runhist.findstep(2:end))/sum(runhist.itersub(2:end))); 
   fprintf('\n scaled-norm(X) = %3.1e, scaled-norm(y) = %3.1e,',...
           normXscaled,normyscaled); 
   fprintf(' scaled-norm(Z) = %3.1e',normZscaled);     
   fprintf('\n norm(X) = %3.1e,        norm(y) = %3.1e,',...
           ops(X,'norm'),norm(y));     
   fprintf('        norm(Z) = %3.1e',ops(Z,'norm')); 
   fprintf('\n--------------------------------------------------------')
   fprintf('------------------------')
   fprintf('\n')
%%
   X(sdpblkidx) = X; Z(sdpblkidx) = Z; 
%%********************************************************************
%%********************************************************************
%% findstep: strong Wolfe line search
%%
%%********************************************************************

   function  [par,y,W,Wp,alp,iter] = ...
              findstep(blk,At,par,y0,W0,Wp0,dy,tol,options);

   if ~exist('tol'); tol = 1e-4; end
   if ~exist('options'); options = 1; end

   maxit = ceil(log(1/(tol+eps))/log(2));
   c1 = 1e-4; c2 = 0.9; 

   b   = par.b; 
   sig = par.sig;
   if (par.use_proximal)
      yhat   = par.yhat;  
      sighat = par.sighat; 
      H      = par.H;  
      H2     = par.H2; 
   end
%%
   g0  = dy'*(b-AXfun(blk,At,Wp0)); 
   Ly0 = b'*y0 - ops(Wp0,'norm')^2/(2*sig);  
   if (par.use_proximal); 
      g0  = g0 - dy'*(H2.*(y0-yhat))/sighat; 
      Ly0 = Ly0 - norm(H.*(y0-yhat))^2/(2*sighat);
   end
   if (par.use_smoothing) 
      Ly0 = Ly0 - sum(par.mu.*par.logphi); 
   end 
   if (g0 <= 0)
      alp = 0; iter = 0; 
      y = y0; W = W0; Wp = Wp0; %%
      fprintf('\n Need an ascent direction, %2.1e  ',g0); 
      return;
   end  
%%
   sigAtdy = Atyfun(blk,At,sig*dy);
   alp = 1; alpconst = 0.5; 
   for iter = 1:maxit
      if (iter==1);          
         alp = 1; LB = 0; UB = 1; 
      else
         alp = alpconst*(LB+UB);
      end
      y = y0 + alp*dy; 
      W = ops(W0,'+',ops(alp,'*',sigAtdy)); 
      [Wp,par] = feval(par.project,blk,W,par);
      normWp = ops(par.phi,'norm'); 
      galp = dy'*(b-AXfun(blk,At,Wp));  
      Ly   = b'*y - normWp^2/(2*sig);
      if (par.use_proximal)
         galp = galp - dy'*(H2.*(y-yhat))/sighat;
         Ly = Ly - norm(H.*(y-yhat))^2/(2*sighat); 
      end   
      if (par.use_smoothing) 
         Ly = Ly - sum(par.mu.*par.logphi); 
      end       
      if (iter==1)
         gLB = g0; gUB = galp; 
         if (sign(gLB)*sign(gUB) > 0)
            fprintf('|')
            return;             
         end
      end
      if (Ly-Ly0-c1*alp*g0 > -1e-10/max(1,abs(Ly))) 
         if (options == 1)
            fprintf(':')
            return;
         elseif (options == 2)
            if (abs(galp) < c2*abs(g0)) & (abs(galp) < tol) 
              fprintf(':')
              return;
            end
         end
      end
      if (sign(galp)*sign(gUB) < 0)
         LB = alp; gLB = galp;
      elseif (sign(galp)*sign(gLB) < 0) 
         UB = alp; gUB = galp; 
      end
   end 
   fprintf('m'); 
%%********************************************************************
