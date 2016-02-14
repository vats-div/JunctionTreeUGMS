%%************************************************************************
%% use the alternating direciton method of multiplier (ADMM) 
%% for the primal problem to generate a good starting point
%%
%%************************************************************************

   function [X,y,Z,par,iter,sig] = ADMMprimal(blk,At,C,b,mu,par,X,y,Z)

   if ~exist('X');   
      X = ops(blk,'identity'); 
      Z = ops(blk,'identity');
      y = zeros(length(b),1);   
   end
   sig = 1;  
   tol = 1e-2; 
   printyes = 1;  
  
   tstart = clock;
   maxiter = 30;
   precond = 1;   
   steplen = 1.618; %% steplength
              
   fprintf('\n *****************************************************');
   fprintf('**********************************');
   fprintf('\n  Alternating direction method of Multiplier (Primal) for generating initial X,y,Z');
   fprintf('\n ***************************************************************');
   fprintf('************************');
   fprintf('\n   it    pinfeas  dinfeas   res_ratio  cg');
   fprintf('   time   sigma');
   fprintf('\n ***************************************************************');
   fprintf('************************');
%%
   sig = 2;
   maxiter = 30; 
   normb = max(1,norm(b)); normC = max(1,ops(C,'norm'));
   V = ops(C,'zeros');
   Y = X;
   Atb = Atyfun(blk,At,b);  
%%   
   for iter = 1:maxiter
      Aty = Atyfun(blk,At,y); 
      tmp = ops(ops(Aty,'+',V),'-',C); %%tmp = ops(Aty,'+',V); 
      rhs = ops(ops(Y,'+',Atb),'+',tmp,1/sig);
      const = 1e-2; 
      tolpsqmr = const*ops(rhs,'norm');
      [X,resnrm,solve_ok] = psqmr_mat(blk,At,rhs,par,tolpsqmr,30);
      par.gamma = mu./sig; 
      W = ops(X,'-',V,1/sig); %%W = ops(X,'-',ops(V,'+',C),1/sig);
      [Y,par] = project3(blk,W,par);
      Rp = b - AXfun(blk,At,Y); 
      step = 1;
      y = y + step*sig*Rp; 
      V = ops(V,'+',ops(Y,'-',X),sig*step);        
      %%
      prim_infeas = norm(Rp)/normb; 
      ttime = etime(clock,tstart); 
      if (rem(iter,5)==1 | iter==maxiter) & (printyes)
         Z = ops(ops(Y,'-',W),'*',sig);
         Aty = Atyfun(blk,At,y);  
         Rd = ops(C,'-',ops(Aty,'+',Z)); 
         dual_infeas = ops(Rd,'norm')/normC; 
         fprintf('\n %3.1d    %3.1e   %3.1e',iter,prim_infeas,dual_infeas); 
         fprintf('   %3.1e   %2.1d ',resnrm(end)/max(1,resnrm(1)),length(resnrm)-1);
         fprintf('   %3.2f   %3.1e',ttime,sig);
      end
      if (max(prim_infeas,dual_infeas) < tol)
         msg = sprintf('prim_infeas < %3.2e',tol); 
         break;
      end
   end
%%*************************************************************************
%%*************************************************************************
%% psqmr: preconditioned symmetric QMR with left (symmetric) preconditioner. 
%%
%% RHS = right-hand side.
%% resnrm = norm of qmr-generated residual vector RHS-linear_operator(X). 
%%
%% linear operator = At*A + I
%%
%% This is for solving the linear system in ADMMprimal.m
%%*************************************************************************

   function  [X,resnrm,solve_ok] = psqmr_mat(blk,At,RHS,par,tol,maxit) 

   rhs = concat(blk,svec(blk,RHS));    
   N = length(rhs); 
   if ~exist('maxit'); maxit = max(50,sqrt(length(rhs))); end;
   if ~exist('printlevel'); printlevel = 1; end;
   if ~exist('tol'); tol = 1e-6*norm(rhs); end;
   if ~exist('L'); par.precond = 0; end
   
   x0 = zeros(N,1);
   solve_ok = 1; 
   precond = 0; 
%%
   x = x0; 
   if (norm(x) > 0) 
      X = smat(blk,unconcat(blk,x));
      AX = AXfun(blk,At,X); 
      tmp = ops(X,'+',Atyfun(blk,At,AX)); 
      Aq = concat(blk,svec(blk,tmp)); 
      r = rhs - Aq;
   else
      r = rhs;
   end
   err = norm(r); resnrm(1) = err; 
%%
   if (precond == 0)
      q = r;
   elseif (precond == 1)      
      q = L.invdiagM.*r;
   end
   tau_old  = norm(q);      
   rho_old  = r'*q; 
   theta_old = 0; 
   d = zeros(N,1); 
   res = r; Ad = zeros(N,1);
%%      
%% main loop
%%
   tiny = 1e-30; 
   for iter = 1:maxit 
       %%
       X = smat(blk,unconcat(blk,q));
       AX = AXfun(blk,At,X); 
       tmp = ops(X,'+',Atyfun(blk,At,AX)); 
       Aq = concat(blk,svec(blk,tmp)); 
       %%
       sigma = q'*Aq; 
       if (abs(sigma) < tiny)
          solve_ok = 2;
          if (printlevel); fprintf('s1'); end
          break;
       else
          alpha = rho_old/sigma; 
          r = r - alpha*Aq;
       end
       if (precond == 0)
          u = r;
       elseif (precond == 1)
          u = L.invdiagM.*r;
       end
       %%
       theta = norm(u)/tau_old; c = 1/sqrt(1+theta^2); 
       tau = tau_old*theta*c;
       gam = (c^2*theta_old^2); eta = (c^2*alpha); 
       d = gam*d + eta*q;
       x = x + d; 
       %%----- stopping conditions ----
       Ad = gam*Ad + eta*Aq;
       res = res - Ad; 
       err = norm(res); resnrm(iter+1) = err; 
       if (err < tol); break; end  
       %%------------------------------
       if (abs(rho_old) < tiny)
          solve_ok = 2;
          if (printlevel); fprintf('s2'); end
          break;
       else
          rho  = r'*u; 
          beta = rho/rho_old; 
          q = u + beta*q; 
       end
       rho_old = rho; 
       tau_old = tau; 
       theta_old = theta; 
   end
   if (iter == maxit); solve_ok = -2; end
%%************************************************************************
%%*************************************************************************
   function xx = concat(blk,x);

   for p = 1:size(x,1)
      idx(p) = length(x{p});
   end
   xx = zeros(sum(idx),1);
   cumidx = [0,cumsum(idx)];
   for p = 1:size(blk,1)
      xx(cumidx(p)+1:cumidx(p+1)) = x{p};
   end
%%*************************************************************************
%%*************************************************************************
   function xx = unconcat(blk,x);

   for p = 1:size(blk,1)
      pblk = blk(p,:);
      if strcmp(pblk{1},'s')
         idx(p) = sum(pblk{2}.*(pblk{2}+1))/2;
      elseif strcmp(pblk{1},'l')
         idx(p) = sum(pblk{2});
      end
   end
   xx = cell(size(blk,1),1);
   cumidx = [0,cumsum(idx)];
   for p = 1:size(blk,1)
      xx{p} = x(cumidx(p)+1:cumidx(p+1));
   end
%%*************************************************************************


