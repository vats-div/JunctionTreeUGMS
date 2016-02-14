%%***************************************************************
%% boundary point method, based on augmented lagrangian 
%% solves: min <C,X> s.t. A*X(:) = b; X psd
%%         max b'y s.t. C - A^t(y) = Z psd
%% maxiter: max number of iterations
%% sigma: penalty parameter for augmented Lagrangian 
%%        a typical value (b and C normalized) is between .1 and 10
%% tol: stopping condition
%%      we stop once primal and dual relative infeasibility < tol
%% 
%%***************************************************************

   function [obj,X,y,Z,runhist,sigma,R,Rt] = mprwmod2(blk,At,C,b,OPTIONS)
 
   if (nargin < 5); OPTIONS = []; end

   tol        = 1e-6;
   maxiter    = 2000; 
   sigma      = 1; 
   scale_data = 0;  
   printlevel = 1; 
   update_it  = 10;  %% update sigma every update_it iterations 
   if isfield(OPTIONS,'tol');        tol        = OPTIONS.tol; end
   if isfield(OPTIONS,'maxiter');    maxiter    = OPTIONS.maxiter; end
   if isfield(OPTIONS,'sigma');      sigma      = OPTIONS.sigma; end 
   if isfield(OPTIONS,'scale_data'); scale_data = OPTIONS.scale_data; end 
   if isfield(OPTIONS,'printlevel'); printlevel = OPTIONS.printlevel; end 
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
      fprintf('\n--------------------------------------------------------')
      fprintf('---------------------');
      fprintf('\n  scale data = %2.0f',scale_data); 
      fprintf('\n--------------------------------------------------------')
      fprintf('---------------------');
      fprintf('\n'); 
   end
%%
%%
%% scaling
%%
   b0 = b; C0 = C; 
   if (scale_data)
      normA0 = min(1e12,max(1,ops(At,'norm'))); 
      normb0 = max(1,norm(b)); normC0 = min(1e12,max(1,ops(C,'norm'))); 
      for p = 1:size(blk,1); At{p} = At{p}/normA0; end
      C = ops(C,'/',normC0); 
      b = b/normb0; 
      objscale = normb0*normC0/normA0; 
   else
      normA0 = 1; normb0 = 1; normC0 = 1;
      objscale = 1;  
   end
   normb = max(1,norm(b)); normC = max(1,ops(C,'norm'));
%%
   tstart = clock; 
   m = length(b); 
   X = ops(C,'zeros'); 
   Z = X;
   y = zeros(m,1);  
%%
%% compute Cholesky factorization of  A*A' 
%%
   [L,AAt] = cholAAt(blk,At,m); 
   if (L.p>0); 
      fprintf(' rows of A linearly dependent. p = %4.0d, m = %4.0d\n',p,m);
      error(' nothing done');
   end
   normX = ops(X,'norm'); 
   Aty = Atyfun(blk,At,y); 
   Rp  = b-AXfun(blk,At,X); 
   Rd  = ops(C,'-',ops(Z,'+',Aty)); 
   normRp = norm(Rp); 
   normRd = ops(Rd,'norm'); 
   obj = objscale*[blktrace(blk,C,X),b'*y]; 
   prim_infeas = normRp/normb; 
   dual_infeas = normRd/normC; 
   runhist.prim_obj(1) = obj(1); 
   runhist.dual_obj(1) = obj(2); 
   runhist.gap(1)      = (obj(1)-obj(2))/(1+sum(abs(obj))); 
   runhist.prim_infeas(1) = prim_infeas;
   runhist.dual_infeas(1) = dual_infeas;
   runhist.cputime(1)     = etime(clock,tstart);
   
   if (printlevel)
      fprintf(' it     time   primal           dual         log(Rp)');
      fprintf(' log(Rd)  sigma   rankX\n');
   end
%%
   for iter = 1:maxiter
      Xold = X; normXold = normX; 
      rhs = AXfun(blk,At,ops(C,'-',Z)) + Rp/sigma;
      y   = linsysolvefun(L,rhs);
      Aty = Atyfun(blk,At,y); 
      M = ops(ops(Aty,'-',C),'+',ops(X,'/',sigma));
      [Mp,Z,d,rankMp] = project(blk,M);
      X = ops(Mp,'*',sigma); 
      normX = ops(X,'norm'); 
      Rp = b-AXfun(blk,At,X); 
      %% Rd = ops(C,'-',ops(Z,'+',Aty)); normRd = ops(Rd,'norm'); 
      normRp = norm(Rp); 
      normRd = sqrt(normX^2 + normXold^2 - 2*blktrace(blk,Xold,X))/sigma;
      obj = objscale*[blktrace(blk,C,X),b'*y]; 
      prim_infeas = normRp/normb; 
      dual_infeas = normRd/normC; 
      runhist.prim_obj(iter+1) = obj(1); 
      runhist.dual_obj(iter+1) = obj(2); 
      runhist.gap(iter+1)      = (obj(1)-obj(2))/(1+sum(abs(obj))); 
      runhist.prim_infeas(iter+1) = prim_infeas;
      runhist.dual_infeas(iter+1) = dual_infeas;
      runhist.cputime(iter+1)     = etime(clock,tstart); 
      if (mod(iter,10)==0) & (printlevel)
         fprintf( '%3.0d %8.1f  %9.8e  %9.8e  %3.2f  %3.2f  %3.2e  %2.0d\n',...
         iter,runhist.cputime(iter+1),obj(1),obj(2),...
         log10(prim_infeas),log10(dual_infeas),sigma,rankMp(1));
      end
      %%
      %% check stopping conditions
      %%
      if (max(prim_infeas,dual_infeas) < tol) | (iter > maxiter)  
         fprintf( '%3.0d %8.1f  %9.8e  %9.8e  %3.2f  %3.2f  %3.2e\n',...
         iter,runhist.cputime(iter+1),obj(1),obj(2),...
         log10(prim_infeas),log10(dual_infeas),sigma);
         if (iter>maxiter)  
            fprintf('max outer iterations reached. \n');
         end
         break; 
      end
      %%
      %% check for reduction of sigma
      %%
      if (prim_infeas > dual_infeas) & (mod(iter,update_it)==0)
         sigma = sigma * 0.9; %% .75
      end
   end
   if (scale_data)
      X = ops(X,'*',normb0/normA0); y = y*(normC0/normA0); 
      Z = ops(Z,'*',normC0);
      for p = 1:size(blk,1); At{p} = At{p}*normA0; end
      Rp = b0 - AXfun(blk,At,X); 
      Rd = ops(C0,'-',ops(Atyfun(blk,At,y),'+',Z));    
      prim_infeas = norm(Rp)/max(1,norm(b0));
      dual_infeas = ops(Rd,'norm')/max(1,ops(C0,'norm'));             
   end
   rel_gap = (obj(1)-obj(2))/max(1,mean(abs(obj)));
   trXZ = blktrace(blk,X,Z); 
   fprintf('\n--------------------------------------------------------')
   fprintf('\n primal objval = %9.8e',obj(1));
   fprintf('\n dual   objval = %9.8e',obj(2));
   fprintf('\n relative gap  = %3.2e',rel_gap);
   fprintf('\n trace(XZ)     = %3.2e',trXZ);
   fprintf('\n prim_infeas   = %3.2e',prim_infeas);
   fprintf('\n dual_infeas   = %3.2e',dual_infeas);
   fprintf('\n CPU time      = %3.1f',runhist.cputime(end)); 
   fprintf('\n norm(X) = %3.1e, norm(y) = %3.1e, norm(Z) = %3.1e',...
           ops(X,'norm'),norm(y),ops(Z,'norm'));     
   fprintf('\n--------------------------------------------------------')
   fprintf('\n')
%%
   X(sdpblkidx) = X; Z(sdpblkidx) = Z; 
%%********************************************************************

