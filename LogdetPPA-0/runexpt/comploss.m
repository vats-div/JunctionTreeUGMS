%%****************************************************************************
%%****************************************************************************

  function [Xsp,Qloss,Zloss,TP,FP,TN,FN] = comploss(Xsp,invSigma); 

  n = length(Xsp); 
  tmp = 2*spones(Xsp)-spones(invSigma);
  TP = length(find(tmp==1))/n^2*100;
  FN = length(find(tmp==-1))/n^2*100; %% TP+FN = nnz(invSigma)/n^2*100
  FP = length(find(tmp==2))/n^2*100;  %% 100-FP-TN = nnz(invSigma)/n^2*100
  TN = length(find(tmp==0))/n^2*100;
  invXsp = inv(Xsp); invXsp = 0.5*(invXsp+invXsp'); 
  Qloss  = norm(invSigma*invXsp-eye(n),'fro')/n; 
  Zloss  = norm(invSigma-Xsp,'fro')/n; 
  fprintf('\n  Qloss = %4.3e, Zloss = %4.3e',Qloss,Zloss) 
%%****************************************************************************
