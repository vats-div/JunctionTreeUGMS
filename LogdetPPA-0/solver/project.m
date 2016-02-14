%%************************************************************************
%% project2: compute projection onto the cone of positive
%%           semidefinite matrices.
%%
%% [Xp,par] = project2(blk,X,par);
%%
%%************************************************************************

    function [Xp,par] = project2(blk,X,par);

    tol = 1e-10;  %% must be small as it will affect gap 
    addtol = 1e-6; 
    Xp = cell(size(X)); 
    par.P1 = cell(size(X));     par.P2 = cell(size(X)); 
    par.Dsch12 = cell(size(X)); 
    par.posidx = cell(size(X)); 
    par.phi = cell(size(X));
    par.phiminus = cell(size(X));
%%
    for p = 1:size(blk,1) 
       pblk = blk(p,:); 
       if strcmp(pblk{1},'s'); 
          blktmp = pblk{1,2}; 
          if (length(blktmp) > 1); 
             %%error(' each cell can only have one block'); 
          end
	  n = sum(blktmp);
          [Xptmp,P,Dsch12,dd,posidx] = projectfun2(X{p},tol);
          %%***** perturbation *****
          Dsch12 = max(addtol,Dsch12); 
          if ~isempty(posidx)
             negidx = setdiff([1:n],posidx); 
             P1 = P(:,posidx); 
             P2 = P(:,negidx);
             phi = dd(posidx); 
             phiminus = abs(dd(negidx)); 
	  else
             P1 = []; 
             P2 = P; 
             phi = []; 
             phiminus = abs(dd); 
          end
       elseif strcmp(pblk{1},'l');
          P1 = []; P2 = []; dd = []; 
          n  = sum(pblk{1,2});
          Xptmp = zeros(n,1); 
          %%***** perturbation *****
          Dsch12 = addtol*ones(n,1); 
          posidx = find(X{p} > tol);  
          if ~isempty(posidx) 
             Xptmp(posidx)  = abs(X{p}(posidx)); 
             Dsch12(posidx) = ones(length(posidx),1);
          end
	  phi = Xptmp;         
          phiminus = abs(Xptmp-X{p}); 
       end
       Xp{p}     = Xptmp;
       par.P1{p} = P1; 
       par.P2{p} = P2; 
       par.posidx{p} = posidx; 
       par.Dsch12{p} = Dsch12; 
       par.phi{p} = phi;
       par.phiminus{p} = phiminus; 
    end
%%***************************************************************************

    function [Xp,V,Dsch12,d,posidx] = projectfun2(X,tol);

    n = length(X);
    if (norm(X-X','fro') > 1e-15*norm(X,'fro'))
       %fprintf(' nonSym '); 
       X = 0.5*(X+X'); 
    end
    if (exist('mexeig')==3)
       [V,D] = mexeig(full(X)); 
    else
       [V,D] = eig(full(X)); 
    end
    d = diag(D);
    [d,idx] = sort(real(d));
    idx = idx(n:-1:1); d = d(n:-1:1); 
    V = V(:,idx); 
    posidx = find(d > tol);
    if isempty(posidx)
       Xp = sparse(n,n); 
       Dsch12 = [];
    elseif (length(posidx) == n)
       Xp = X; 
       Dsch12 = []; 
    else
       r = length(posidx); s = n-r; 
       negidx = [r+1:n];      
       dp = abs(d(posidx)); 
       dn = abs(d(negidx));
       Vtmp = V(:,posidx)*diag(sqrt(dp)); 
       Xp = Vtmp*Vtmp'; 
       Xp = 0.5*(Xp+Xp'); 
       Dsch12 = (dp*ones(1,s))./(dp*ones(1,s) + ones(r,1)*dn');  
    end
%%************************************************************************
