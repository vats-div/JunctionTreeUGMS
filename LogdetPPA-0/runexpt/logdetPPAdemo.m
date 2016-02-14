%%**********************************************************************
%% This subscript is to test the codes implementing 
%% a primal-proximal point method for log-determinant 
%% optimization problem
%%
%% min   <Sigma, X> -logdet(X) + sum_{ij} rho_{ij} |X_{ij}|
%% s.t.  X_{ij} = 0, (i,j) in Omega,
%%       X is positive definite,
%% where Omega := {(I(k), J(k)): 1 <= k <= |I|}
%%
%% Such a problem arises from covariance selection problem in statistics.
%% Detials on the numerical experiments conducted for such problems can be found
%% in the paper given in Citation.txt. 
%%**********************************************************************
%% LogdetPPA: version 0
%% Copyright (c) 2009 by
%% Chengjing Wang, Defeng Sun, and Kim-Chuan Toh
%% Last modified: 27 May 1010
%%**********************************************************************

  warning off
  clear all
  HOME = './LogdetPPA'; 
  addpath(strcat(HOME,'/solver/'))
  addpath(strcat(HOME,'/solver/mexfun'))
  addpath(strcat(HOME,'/util/'))
%%
  fprintf('\n----------------------------------------------------------');
  fprintf('\n An example to demonstrate how to use LogdetPPA to solve a'); 
  fprintf('\n log-determinant SDP arising from covariance selection'); 
  fprintf('\n problem in statistics.\n'); 
%%
%%
%%
  ZERO     = 1e-3;  %% small perturbation
  signoise = 0.15;  %% Signal to noise ratio
  density  = 0.1;
%%
%% L_options = 0 without L1 term in the objective
%%           = 1 with L1 term in the objective
%%
  L_options = 0; 
  invSigma_options = 1; %% default = 1; 
  S_options        = 2; %% default = 2; 
  dim = [100,200,500];  
%%
  for dim_index = [1];  
     rand('state',0); 
     randn('state',0);

     ttime = clock;
     n = dim(dim_index);
     sstr = ['-L',num2str(L_options),'-S',num2str(S_options)]; 
     if (invSigma_options==1)
        fnamestr = [num2str(n),sstr]; 
     else
        fnamestr = [num2str(n),sstr]; 
        fnamestr = [fnamestr,'-invS',num2str(invSigma_options)]; 
     end
     if (invSigma_options==1) %% used by d'Aspremont
        numnz = max(10,ceil(density*n^1.68)); 
        len = max(10,numnz); 
        ii = ceil(n*rand(len,1)); 
        jj = ceil(n*rand(len,1)); 
        vv = sign(rand(len,1)-0.5); 
        A = spconvert([ii,jj,vv;n,n,0]); 
        A = A'*A; A = 0.5*(A+A'); 
        diagA = diag(A); 
        A2 = A-spdiags(diagA,0,n,n); 
        A2 = max(min(A2,1),-1); 
        A = A2 + spdiags(1+diagA,0,n,n); 
        rho = 1/n^1.5; 
        if (n==30); rho = 0.03; end
     elseif (invSigma_options==2)
        A = sprandsym(n,density);
        rho = 1/n^1.5;
     end
     A = full(A); 
     [eigvec,eigval] = mexeig(A); 
     mineigA  = min(diag(eigval));
     invSigma = A + (max(-1.2*mineigA,0)+ZERO)*eye(n); 
     [eigvec,eigval] = mexeig(invSigma); eiginvSigma = full(diag(eigval)); 
     fprintf('\n nnz(invSigma)/n^2 = %3.2e',nnz(invSigma)/n^2*100); 
     fprintf('\n max(eiginvSigma), min(eiginvSigma) = %3.2e, %3.2e',...
              max(eiginvSigma),min(eiginvSigma)); 
     %%
     %% Iomega, Jomega
     %% 
     cutoff = 5; 
     [II,JJ] = find(invSigma==0); 
     idx = find(JJ-II >= cutoff); 
     Iomega = II(idx); 
     Jomega = JJ(idx);
     %%
     %% sample covariance matrix S
     %%
     if (S_options==1) %% used by d'Aspremont
        Pert  = 2*rand(n,n)-1;  
        Pert  = 0.5*(Pert+Pert');
        normratio = norm(inv(invSigma),'fro')/norm(Pert,'fro'); 
        %%normratio = 1;
        S = inv(invSigma) + signoise*normratio*Pert; 
        S = 0.5*(S+S');  
        mineigS = min(eig(S));       
        S = S + (max(-mineigS,0)+ZERO)*eye(n); 
        Qloss_LW = inf;
        Zloss_LW = inf;
     elseif (S_options==2) %% statistically more realistic
        [Q,D] = mexeig(full(invSigma)); 
        eigSigma = 1./max(1e-8,full(diag(D)));  
        Q = Q*spdiags(sqrt(eigSigma),0,n,n);
        ss = ceil(2*n); %% sample size
        Samp = zeros(n,ss); 
        for k = 1:ss
           Samp(:,k) = Q*randn(n,1); 
        end
        S = cov(Samp');  
        trS = sum(diag(S)); 
        beta1 = (ss-2)/ss*norm(S,'fro')^2 + trS^2; 
        beta2 = (ss+2)*(norm(S,'fro')^2 - trS^2/n); 
        beta  = min(1,beta1/beta2); 
        LW    = (beta*(trS/n))*eye(n) + (1-beta)*S; 
        Qloss_LW = norm(invSigma*LW-eye(n),'fro')/n; 
        Zloss_LW = norm(invSigma-inv(LW),'fro')/n; 
     end
     [eigvec,eigval] = mexeig(full(S)); 
     eigS = full(diag(eigval)); 
     %%
     m = length(Iomega); 
     fprintf('\n set up data time = %3.2f',etime(clock,ttime)); 
     fprintf('\n length(Iomega) = %2.0d',m);
     fprintf('\n----------------------------------------------------------\n');
     pause(3); 
%%
     [blk,At,C,b,mu] = genSDPdata(Iomega,Jomega,S,L_options,rho); 
     if (L_options==1)
        OPTIONS.sigma      = 10;
        OPTIONS.scale_data = 2; 
     else
        OPTIONS.sigma      = 1; 
        OPTIONS.scale_data = 0;
     end
     OPTIONS.plotyes = 0; 
     OPTIONS.tol     = 1e-6;
     OPTIONS.smoothing = 1;
     [obj,X,y,Z,info,runhist] = logdetPPA(blk,At,C,b,mu,OPTIONS);
     X1 = X{1}; 
     [ii,jj] = find(abs(X1) < 1e-6*max(max(abs(X1))));     
     Xsp = X1.*(1-spconvert([ii,jj,ones(length(ii),1);n,n,0]));
     [Xsp,Qloss,Zloss,TP,FP,TN,FN] = comploss(Xsp,invSigma); 
     fprintf('\n  Qloss_LW = %4.3e, Zloss_LW = %4.3e\n',Qloss_LW,Zloss_LW); 
     runhist.rho              = rho; 
     runhist.signoise         = signoise;
     runhist.nnzinvSigma      = nnz(invSigma); 
     runhist.S_options        = S_options; 
     runhist.invSigma_options = invSigma_options; 
     runhist.eigS             = eigS; 
     runhist.eiginvSigma      = eiginvSigma; 
     runhist.TP = TP; runhist.FN = FN; runhist.FP = FP; runhist.TN = TN; 
     runhist.Qloss = Qloss; runhist.Zloss = Zloss; 
     runhist.Qloss_LW = Qloss_LW; runhist.Zloss_LW = Zloss_LW; 
  end  
%%**********************************************************************  

