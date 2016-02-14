%%**********************************************************************
%% generate SDP data for the following problem: 
%%
%% L_options = 1:
%% min <S,X> - log(det(X)) + rho*sum_{(ij)\not\in Omega} |Xij|
%% s.t   Xij = 0, (ij)\in Omega. 
%% 
%% L_options = 0: 
%% min <S,X> - log(det(X)) 
%% s.t   Xij = 0, (ij)\in Omega. 
%%**********************************************************************

     function [blk,At,C,b,mu] = genSDPdata(Iomega,Jomega,S,L_options,rho); 

     n = length(S); 
     C{1} = S; 
     blk{1,1}= 's'; blk{1,2} = n;
     m = length(Iomega);
     b = zeros(m,1);
     n2 = n*(n+1)/2; 
     %%
     %% generate At corresponding to X_Omega  
     %%
     if (length(Iomega))
        Itmp = Iomega + Jomega.*(Jomega-1)/2; 
        At{1} = spconvert([Itmp,[1:m]',ones(m,1); n2,m,0]); 
     else
	At{1} = [];
     end
     %%
     if (L_options==1)
         [Iall,Jall] = find(triu(ones(n)));
         if (length(Iomega))
             tmp = setdiff([Iall,Jall],[Iomega,Jomega],'rows');
         else
             tmp = [Iall,Jall];
         end
         m2 = size(tmp,1);
         Icomp = tmp(:,1); Jcomp = tmp(:,2);
         Itmp = Icomp + Jcomp.*(Jcomp-1)/2;
         Atmp = spconvert([Itmp,[1:m2]',ones(m2,1); n2,m2,0]);
         blk{2,1} = 'l'; blk{2,2} = 2*m2;
         At{1} = [At{1},Atmp];
         Identity = speye(m2);
         At{2,1} = [sparse(m,2*m2); -Identity,Identity]';
         b = [b; zeros(m2,1)];
         idx = find(Icomp == Jcomp);
         ee = sqrt(2)*ones(m2,1);
         if ~isempty(idx); ee(idx) = ones(length(idx),1); end
         C{2,1} = rho*[ee; ee];
         mu = [1; 1e-8];
     elseif (L_options==0)
         mu = 1;
     end
