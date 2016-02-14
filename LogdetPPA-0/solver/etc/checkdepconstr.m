%%*****************************************************************************
%% checkdepconst: compute AAt to determine if the 
%%             constraint matrices Ak are linearly independent. 
%%              
%% [At,b,y,idxB,neardepconstr,feasible,AAt] = checkdepconstr(blk,At,b,y,rmdepconstr);
%%
%% rmdepconstr = 1, if want to remove dependent columns in At
%%             = 0, otherwise.
%% 
%% idxB = indices of linearly independent columns of original At.
%% neardepconstr = 1 if there is nearly dependent columns in At
%%            = 0, otherwise.
%% Note: the definition of "nearly dependent" is dependent on the 
%%       threshold used to determine the small diagonal elements in 
%%       the LDLt factorization of A*At. 
%% 
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*****************************************************************************

   function  [At,b,idxB,neardepconstr,feasible,AAt] = ...
              checkdepconstr(blk,At,b,rmdepconstr);
   
   printlevel   = 1; 
   existlowrank = 0; 
   matlabversion = sscanf(version,'%f');
   matlabversion = matlabversion(1);
   if strcmp(computer,'PCWIN64') | strcmp(computer,'GLNXA64')
      computer_model = 64; 
   else
      computer_model = 32; 
   end
%%  
%% compute AAt
%%
   m = length(b);  
   AAt = sparse(m,m);  
   for p = 1:size(blk,1)
      AAt = AAt + At{p,1}'*At{p,1}; 
   end
   pertdiag = 1e-13*ones(m,1); 
   AAt = AAt + spdiags(pertdiag,0,m,m);
%%
   if (nnz(AAt) < 0.2*m*m); use_spchol=1; else; use_spchol=0; end
   if (use_spchol)
      L.matfct_options = 'spcholmatlab'; 
      [L.R,L.p,L.perm] = chol(AAt,'vector'); 
      L.Rt = L.R';       
      indef = L.p;
   else
      if issparse(AAt); AAt = full(AAt); end;           
      L.matfct_options = 'chol';
      L.perm = [1:m]; 
      [L.R,indef] = chol(AAt); 
   end
   L.d = diag(L.R).^2; 
%%
%% 
%%
   feasible = 1; neardepconstr = 0; 
   if ~issparse(AAt); AAt = sparse(AAt); end 
   if (exist('mexnnz') == 3)
      nnzmatold = mexnnz(AAt);
   else
      nnzmatold = nnz(AAt);
   end
   if (indef) 
      msg = 'AAt is not pos. def.'; 
      idxB = [1:m]; 
      neardepconstr = 1; 
      if (printlevel); fprintf('\n checkdepconstr: %s',msg); end
      return; 
   end
%%
%% find independent rows of A
%%
   dd = zeros(m,1); 
   idxB = [1:m]';
   dd(L.perm) = abs(L.d); 
   idxN = find(dd < 1e-10*mean(L.d));
   ddB = dd(setdiff([1:m],idxN));
   ddN = dd(idxN);
   if ~isempty(ddN) & ~isempty(ddB) & (min(ddB)/max(ddN) < 10) 
      %% no clear separation of elements in dd
      %% do not label constraints as dependent
      idxN = []; 
   end
numdencol = 0; 
   if ~isempty(idxN)     
      neardepconstr = 1; 
      if (printlevel)
         fprintf('\n number of nearly dependent constraints = %1.0d',length(idxN)); 
      end
      if (numdencol==0)
         if (rmdepconstr)
            idxB = setdiff([1:m]',idxN);
            if (printlevel)
               fprintf('\n checkdepconstr: removing dependent constraints...');
            end
            [W,resnorm] = findcoeffsub(blk,At,idxB,idxN);
   	    tol = 1e-8;
            if (resnorm > sqrt(tol))
               idxB = [1:m]'; 
               neardepconstr = 0; 
               if (printlevel)
                  fprintf('\n checkdepconstr: basis rows cannot be reliably identified,'); 
                  fprintf(' abort removing nearly dependent constraints'); 
               end
               return; 
            end
            tmp = W'*b(idxB) - b(idxN);
            nnorm = norm(tmp)/max(1,norm(b)); 
            if (nnorm > tol) 
               feasible = 0; 
               if (printlevel)
                  fprintf('\n checkdepconstr: inconsistent constraints exist,');
                  fprintf(' problem is infeasible.');
               end
            else
               feasible = 1; 
               for p = 1:size(blk,1) 
                  At{p,1} = At{p,1}(:,idxB);
               end
               b = b(idxB);
               AAt = AAt(idxB,idxB);               
            end
	 else
            if (printlevel)
               fprintf('\n there are dependent constraints,');
            end
         end
      else
         if (printlevel)
            fprintf('\n warning: the sparse part of AAt may be nearly singular.');
         end
      end
   end
%%*****************************************************************************
%% findcoeffsub: 
%%
%% [W,resnorm] = findcoeffsub(blk,At,idXB,idXN);
%% 
%% idXB = indices of independent columns of At. 
%% idxN = indices of   dependent columns of At.
%% 
%% AB = At(:,idxB); AN = At(:,idxN) = AB*W
%%
%% SDPT3: version 3.0
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%*****************************************************************************

   function [W,resnorm] = findcoeffsub(blk,At,idxB,idxN);

   AB = []; AN = [];
   for p = 1:size(blk,1) 
      AB = [AB; At{p,1}(:,idxB)];
      AN = [AN; At{p,1}(:,idxN)];
   end
   [m,n] = size(AB); 
%%
%%-----------------------------------------
%% find W so that AN = AB*W
%%-----------------------------------------
%% 
   [L,U,P,Q] = lu(sparse(AB));    
   rhs  = P*AN;
   Lhat = L(1:n,:); 
   W = Q*( U \ (Lhat \ rhs(1:n,:))); 
   resnorm = norm(AN-AB*W,'fro')/max(1,norm(AN,'fro'));
%%*****************************************************************************
