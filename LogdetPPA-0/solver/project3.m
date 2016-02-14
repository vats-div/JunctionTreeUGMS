%%***************************************************************************
%% project3: compute projection onto the cone of positive
%%           semidefinite matrices, with smoothing. 
%%
%% Xp = argmin {0.5*||Y-X||^2 - gamma*logdet(Y) : Y pd}
%% 
%% [Xp,par] = project3(blk,X,par);
%%
%%***************************************************************************

    function [Xp,par] = project3(blk,X,par);

    tol = 1e-10; 
    addtol = 1e-6; 
    gamma = par.gamma; 
    if (length(gamma)==1); gamma = gamma*ones(length(X),1); end

    Xp = cell(size(X)); 
    par.P = cell(size(X));   
    par.posidx = cell(size(X)); 
    par.Dsch = cell(size(X)); 
    par.phi = cell(size(X));
    par.phiminus = cell(size(X));
%%
    for p = 1:size(blk,1) 
       pblk = blk(p,:); 
       if strcmp(pblk{1},'s'); 
          blktmp = pblk{1,2}; 
          if (length(blktmp) > 1); 
             %% error(' each cell can only have one block'); 
             %% treat multiple sub-blocks as a single large block
          end
          n = length(X{p});
          if (norm(X{p}-X{p}','fro') > 1e-15)
             X{p} = 0.5*(X{p}+X{p}'); 
          end
	  if exist('mexeig')
             [P,D] = mexeig(full(X{p}));
          else
	     [P,D] = eig(full(X{p})); 
          end
          d = diag(D); 
          [dd,idx] = sort(real(d));
          idx = idx(n:-1:1); dd = dd(n:-1:1); 
          P = P(:,idx); 
          posidx = find(dd > tol);
          ss = sqrt(dd.*dd + 4*gamma(p)); 
          phi = 0.5*(dd + ss) + eps; 
	  phiminus = max(0.5*(ss-dd),0) + eps; 
          Xptmp = (P*spdiags(phi,0,n,n))*P'; 
          Xptmp = 0.5*(Xptmp+Xptmp'); 
          ee = ones(n,1); 
          Dsch = 0.5*(1 + (dd*ee'+ee*dd')./(ss*ee' + ee*ss'));  
       elseif strcmp(pblk{1},'l');
          P = []; dd = []; 
          posidx = find(X{p} > tol);  
          ss = sqrt(X{p}.*X{p} + 4*gamma(p)); 
          Xptmp = 0.5*(X{p} + ss) + eps;
          phi   = Xptmp;  
          phiminus = max(0.5*(ss-X{p}),0) + eps; 
          Dsch = 0.5*(1 + X{p}./ss); 
       end
       Dsch = Dsch + addtol; 
       Xp{p} = Xptmp; 
       par.P{p} = P;
       par.Dsch{p} = Dsch;
       par.logphi(p,1) = sum(log(phi)); 
       par.phi{p}      = phi; 
       par.phiminus{p} = phiminus; 
       par.posidx{p} = posidx;
    end
%%***************************************************************************


