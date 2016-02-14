%%*************************************************************************
%% psqmr_newton:  symmetric QMR. 
%%
%% b = rhs vector.
%% resnrm = norm of qmr-generated residual vector b-Ax. 
%%*************************************************************************

   function  [dX,dy,resnrm,num_inneriter,solve_ok] = ...
              psqmr_newton(matvecfname,blk,At,borg,par,L,tol,maxit,AL) 
  
   if iscell(borg); 
      b = concat(blk,svec(blk,borg)); 
   else
      b = svec(blk,borg); 
   end
   N = length(b); 
   if ~exist('maxit'); maxit = max(50,sqrt(length(b))); end;
   if ~exist('tol'); tol = 1e-6*norm(b); end; 
   if ~exist('AL'); AL = []; end

   solve_ok = 1; 
%%
   if isfield(par,'P');  par.Pt = ops(par.P,'transpose'); end
   if isfield(par,'P1'); par.P1t = ops(par.P1,'transpose'); end
   if isfield(par,'P2'); par.P2t = ops(par.P2,'transpose'); end
   if isfield(par,'Dsch12'); par.Dsch21 = ops(par.Dsch12,'transpose'); end
%%
   x = zeros(N,1); 
   r = b;  
   normr = norm(r); 
   err = normr; resnrm(1) = err; 
%%
   q = r; 
   tau_old  = normr;      
   rho_old  = normr^2; 
   theta_old = 0; 
   d = zeros(N,1); 
   res = r; Ad = zeros(N,1);
%%      
%% main loop
%%
   tiny = 1e-30; 
   for iter = 1:maxit 
  
       H = smat(blk,unconcat(blk,q)); 
       Y = DPhi(blk,par,H);
       rhs = -AXfun(blk,At,Y); 
       tol_dy = 5e-3*norm(rhs);  %% old=1e-2; 
       [dy,resnrmtmp] = psqmr(matvecfname,blk,At,rhs,par,L,tol_dy,10*max(2,iter));
       num_inneriter(iter) = length(resnrmtmp)-1; 
       Y = ops(H,'+',Atyfun(blk,At,par.sig*dy)); 
       Aq = concat(blk,svec(blk,ops(H,'-',DPhi(blk,par,Y)))); 
       %%
       sigma = q'*Aq; 
       if (abs(sigma) < tiny)
          solve_ok = 2; 
          %fprintf('s1'); 
          break;
       else
          alpha = rho_old/sigma; 
          r = r - alpha*Aq;
       end
       normr = norm(r);  
       %%
       theta = normr/tau_old; c = 1/sqrt(1+theta^2); 
       tau = tau_old*theta*c;
       gam = (c^2*theta_old^2); eta = (c^2*alpha); 
       d = gam*d + eta*q;
       x = x + d; 
       %%----- stopping conditions ----
       Ad = gam*Ad + eta*Aq;
       res = res - Ad; 
       err = norm(res); resnrm(iter+1) = err; 
       if (err < tol); break; end  
       %%----------------------------- 
       if (abs(rho_old) < tiny)
          solve_ok = 2; 
          %fprintf('s2');
          break;
       else
          rho  = normr^2; 
          beta = rho/rho_old; 
          q = r + beta*q; 
       end
       rho_old = rho; 
       tau_old = tau; 
       theta_old = theta; 
   end
   if (solve_ok ~= -1); %fprintf(' '); 
   end
   if iscell(borg); xtmp = unconcat(blk,x); clear x; x = xtmp; end
   dX = smat(blk,x); 
%%*************************************************************************
%% compute P*(Dsch.*(Pt*H*P))*Pt
%%*************************************************************************
   function Y = DPhi(blk,par,H); 
 
   Y = cell(size(blk,1),1); 
   if isfield(par,'P')
      for p = 1:size(blk,1);
         pblk = blk(p,:);
         if strcmp(pblk{1},'s')
            tmp = par.Pt{p}*H{p}*par.P{p}; 
            tmp = par.Dsch{p}.*tmp; 
            tmp = par.P{p}*tmp*par.Pt{p}; 
            Y{p} = 0.5*(tmp+tmp');
         elseif strcmp(pblk{1},'l')
            Y{p} = par.Dsch{p}.*H{p}; 
         end
      end
   elseif isfield(par,'P1')
      for p = 1:size(blk,1);
         pblk = blk(p,:);
         n = sum(pblk{2}); 
         if strcmp(pblk{1},'s')
            Y{p} = sparse(n,n); 
            rr = size(par.P1{p},2);
            if (rr > 0 & rr < n) 
               if (rr <= n/2) 
                  tmp0 = par.P1t{p}*H{p};
                  tmp1 = (tmp0*par.P1{p})*par.P1t{p};         
                  tmp2 = par.Dsch12{p}.*(tmp0*par.P2{p});
                  tmp2 = tmp2*par.P2t{p}; 
                  tmp3 = par.P1{p}*(0.5*tmp1 + tmp2);
                  Y{p} = tmp3+tmp3';
               else
                  tmp0 = par.P2t{p}*H{p};
                  tmp1 = (tmp0*par.P2{p})*par.P2t{p};         
                  tmp2 = (1-par.Dsch21{p}).*(tmp0*par.P1{p});
                  tmp2 = tmp2*par.P1t{p}; 
                  tmp3 = par.P2{p}*(0.5*tmp1 + tmp2);
                  Y{p} = H{p}-tmp3-tmp3'; 
               end
            elseif (rr == n)
               Y{p} = H{p};      
            end
         elseif strcmp(pblk{1},'l')
            Y{p} = sparse(n,1); 
            if (~isempty(par.Dsch12{p}))
               Y{p} = par.Dsch12{p}.*H{p}; 
            end
         end
      end
   end   
%%*************************************************************************
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
